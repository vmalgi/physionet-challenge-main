'''This example demonstrates the use of Convolution1D for text classification.

Run on GPU: THEANO_FLAGS=mode=FAST_RUN,device=gpu,floatX=float32 python imdb_cnn.py

Get to 0.835 test accuracy after 2 epochs. 100s/epoch on K520 GPU.
'''

from __future__ import print_function
import numpy as np
np.random.seed(1337)  # for reproducibility
from scipy import stats

from keras.preprocessing import sequence
from keras.models import Sequential
from keras.optimizers import SGD
from keras.layers.core import Dense, Dropout, Activation, Flatten
from keras.layers.embeddings import Embedding
from keras.layers.convolutional import Convolution2D, MaxPooling2D
#from keras.layers.normalization import BatchNormalization
from keras.utils import np_utils
from keras.regularizers import l2, activity_l2
from keras.layers.normalization import BatchNormalization
from keras.layers.advanced_activations import ELU
from keras.callbacks import ModelCheckpoint
from keras.constraints import maxnorm

import hbdata_predict as dataimport
from keras.models import model_from_json
import sys


# set parameters:
batch_size = 64

def write_answer(filename, result, resultfile="answers.txt"):
        fo = open(resultfile, 'a')
        fo.write(str(filename) + "," + str(result) + "\n")
        fo.close()

        return True


print('Build model...')
model = model_from_json(open('hb_model_orthogonal_experiment_norm.json').read())
model.load_weights('hb_weights_orthogonal_experiment_norm.hdf5')
sgd = SGD(lr=0.00001, decay=1e-6, momentum=0.9)

model.compile(loss='categorical_crossentropy', optimizer=sgd)

filename = sys.argv[1]

(X, labels) = dataimport.data_from_file(filename=str(filename)+".wav", label=0, width=129, height=256, max_frames=10)
predictions = np.zeros(len(X))
z = 0
for frame in X:
	predict_frame = np.zeros((1, 3, 129, 129))
	predict_frame[0] = frame
	predictions_all = model.predict_proba(predict_frame, batch_size=batch_size)
	predictions[z] = predictions_all[0][1]
	z += 1
average = np.average(predictions)
average_prediction = round(average)
		
if int(average_prediction) == 0.0:
		#append file with -1
	write_answer(filename=filename, result="-1")
else:
		#append file with 1
	write_answer(filename=filename, result="1")
	
