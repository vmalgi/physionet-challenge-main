#!/usr/bin/python -tt

# # PhysioNet Test Entry
# # This script will be invoked by next.sh
# ### Input: a PCG record to process (.wav). 
# ### Output: appends classification result to answers.txt

import sys
import numpy as np
import pandas as pd
import sklearn
import tensorflow as tf

from sklearn import preprocessing
import datetime
import scipy.io as sio
import pickle
import math
import os

import scipy.signal as sig
from features import mfcc
import scipy.io.wavfile as wav

# Stay within time limit
MAX_SEGMENTS = 200

num_features = 6
duration = 3
window_step = 0.01
Fs = 2000.0

do_feature_scaling = True

frame_length = int(duration / window_step)

# Load and start Matlab Engine

import matlab.engine
eng = matlab.engine.start_matlab()
eng.addpath(r'matlab/', nargout=0)

# Saved Model - Convolutional Neural Network - Mel Frequency Cepstral Coefficients
saved_model = 'cnn'

results_loc = '.'
    
# Deserialize scaler
scaler = np.load('saved_models/{}/scaler.npy'.format(saved_model))
    
# Convnet Parameters

num_channels = 1
num_labels = 2

patch_size1 = 20
patch_size2 = 10

patch_height1 = 2
patch_height2 = 2

max_pool_length1 = 20
max_pool_length2 = 4

max_pool_stride1 = 5
max_pool_stride2 = 2

depth = 64
num_hidden1 = 1024
num_hidden2 = 512

## Reconstruct Graph
   
graph = tf.Graph()

random_seed = 5001
tf.set_random_seed(random_seed)

with graph.as_default():
    
    with tf.device('/gpu:5'):
        
        # Input data.
        with tf.name_scope('initial_data'):
            # Add test dataset
            tf_test_dataset = tf.placeholder(tf.float32, shape=(None, num_features, frame_length, num_channels))        
            tf_batch_size = tf.placeholder(tf.int32)
        
        def reshape(original_tensor, batch_size):
            return tf.reshape(original_tensor, tf.pack([batch_size, -1]))
        
        # Variables.

        # Dimensions for conv weights are: 
        # patch_height x patch_width x #channels x depth
        with tf.name_scope('layer1'):
            layer1_weights = tf.Variable(tf.truncated_normal([patch_height1, patch_size1, num_channels, depth], stddev=0.1))
            layer1_biases = tf.Variable(tf.zeros([depth]))

        with tf.name_scope('layer2'):
            layer2_weights = tf.Variable(tf.truncated_normal([patch_height2, patch_size2, depth, depth], stddev=0.1))
            layer2_biases = tf.Variable(tf.constant(1.0, shape=[depth]))

        # Unroll convnet into feature vector (layer_size * depth)
        with tf.name_scope('fully_connected_layer'):
            layer3_weights = tf.Variable(tf.truncated_normal([(frame_length // max_pool_stride1 // max_pool_stride2) * num_features * depth, num_hidden1], stddev=0.1))
            layer3_biases = tf.Variable(tf.constant(1.0, shape=[num_hidden1]))

            layer4_weights = tf.Variable(tf.truncated_normal([num_hidden1, num_hidden2], stddev=0.1))
            layer4_biases = tf.Variable(tf.constant(1.0, shape=[num_hidden2]))

            layer5_weights = tf.Variable(tf.truncated_normal([num_hidden2, num_labels], stddev=0.1))
            layer5_biases = tf.Variable(tf.constant(1.0, shape=[num_labels]))

        # Model.
        def model(data):
            # Dimensions for strides are:
            # batch x patch_height x patch_width x #channels
            # e.g. [1, 2, 2, 1]

            print data.get_shape().as_list()

            with tf.name_scope('conv1'):
                conv = tf.nn.conv2d(data, layer1_weights, [1, 1, 1, 1], padding='SAME')
                bias = tf.nn.bias_add(conv, layer1_biases)
                conv1 = tf.nn.relu(bias)

            # Do max-pooling with patch size of 1x2 and stride of 1x2, include all batches and channels
            with tf.name_scope('max_pool1'):
                maxpool1 = tf.nn.max_pool(conv1, [1, 1, max_pool_length1, 1], [1, 1, max_pool_stride1, 1], padding='SAME')

            print maxpool1.get_shape().as_list()

            with tf.name_scope('conv2'):
                conv = tf.nn.conv2d(maxpool1, layer2_weights, [1, 1, 1, 1], padding='SAME')
                bias = tf.nn.bias_add(conv, layer2_biases)
                conv2 = tf.nn.relu(bias)
                
            with tf.name_scope('max_pool2'):
                maxpool2 = tf.nn.max_pool(conv2, [1, 1, max_pool_length2, 1], [1, 1, max_pool_stride2, 1], padding='SAME')

            shape = maxpool2.get_shape().as_list()
            print shape

            # Fully connected layer, batch_size x total_features
            # Rollout height, width and feature_map into total features
            reshaped = reshape(maxpool2, tf_batch_size)

            with tf.name_scope('mlp_hidden'):
                hidden = tf.nn.relu(tf.matmul(reshaped, layer3_weights) + layer3_biases)
                hidden = tf.nn.relu(tf.matmul(hidden, layer4_weights) + layer4_biases)

            with tf.name_scope('mlp_y'):
                y = tf.matmul(hidden, layer5_weights) + layer5_biases
                
            return y
                
        # Predictions for the test dataset
        test_prediction = tf.nn.softmax(model(tf_test_dataset))
        saver = tf.train.Saver()     
                
# Reformat data

def reformat(dataset, labels):
    dataset = dataset.reshape((-1, num_features, frame_length, num_channels)).astype(np.float32)
    labels = (np.array([-1, 1]) == labels[:,None]).astype(np.float32)
    return dataset, labels
    
def accuracy(predictions, labels):
    return (100.0 * np.sum(np.argmax(predictions, 1) == np.argmax(labels, 1)) / predictions.shape[0])
    
     
def make_prediction(filename):
	with tf.Session(graph=graph, config=tf.ConfigProto(allow_soft_placement=True, log_device_placement=False)) as session:

		tf.initialize_all_variables().run()

		# Restore and initialize the saved model
		ckpt = tf.train.get_checkpoint_state('saved_models/{}'.format(saved_model))
		if ckpt and ckpt.model_checkpoint_path:

			saver.restore(session, ckpt.model_checkpoint_path)
			print('Checkpoint found - model restored')

			# Evaluate instance
			# -----------------

			# Call matlab routine for segmenting PCG data
			A, X = eng.generateSegmentedDataset('./', filename, nargout=2)

			# Convert back to numpy format
			A = np.array(A)

			# Get start of heart beat indices (S1) in seconds
			# Matlab resampled to 1KHz, so multiply by 2 (assumes 2Khz sampling rate)
			s1_indices = A[:,0] * 2 / Fs
			frame_starts = s1_indices / window_step

			# Derive mfcc features for entire recording
			(rate, sig) = wav.read('{}.wav'.format(filename))
			mfcc_feat = mfcc(sig, rate, winstep=window_step, numcep=num_features)

			X_test = []

			curr_frame = -1
			count = 0
			for frame_start in frame_starts:
				if count < MAX_SEGMENTS:
					curr_frame = round(frame_start)
					mfcc_slice = mfcc_feat.T[:,int(curr_frame):int(curr_frame+frame_length)]

					# Only include full duration slices
					if mfcc_slice.shape[1] == frame_length:
						X_test.append(mfcc_slice)
				count += 1

			# Ensure at least one 5 second window exists
			if len(X_test) > 0:

				if do_feature_scaling:

					# Reshape data for scaling
					test_dataset = []
					for x in X_test:
						feats = []
						for i in range(num_features):
							feats.append(x[i, :])
						test_dataset.append(np.hstack(feats).squeeze())
					test_dataset = np.vstack(test_dataset)
					
					# Apply scaler
					X_test =  (test_dataset - scaler[0]) / scaler[1]

				test_dataset, test_labels = reformat(X_test, np.ones(len(X_test)))   # Dummy labels 
				print('Test set', test_dataset.shape, test_labels.shape)

				# Predict with loaded model

				# Feed to loaded model
				feed_dict = {tf_test_dataset : test_dataset, tf_batch_size : len(X_test)}
				predictions = session.run(test_prediction, feed_dict=feed_dict)

				# Ensemble strategy: maximize average
				avg_pred = np.average(predictions, axis=0)
				pred = 1 if np.argmax(avg_pred) == 1 else -1

				# Write to predictions.txt
				#f = open('{}/predictions.txt'.format(results_loc), 'a')
				#f.write('{},{},{}\n'.format(filename,avg_pred[0],avg_pred[1]))
				#f.close()

			else:
				pred = 0   # When all else fails, predict uncertain
				#f = open('{}/predictions.txt'.format(results_loc), 'a')
				#f.write('{},0.5, 0.5\n'.format(filename))
				#f.close()

			# Write to answer.txt
			f = open('{}/answers.txt'.format(results_loc), 'a')
			f.write('{},{}\n'.format(filename,pred))
			f.close()
		else:
			print 'Model could NOT be restored'
	
                
def main():
	if len(sys.argv) != 2:
		print 'usage: ./predict.py filename'
		sys.exit(1)
	
	filename = sys.argv[1]
	make_prediction(filename)
	eng.quit()
	
if __name__ == '__main__':
	main()


