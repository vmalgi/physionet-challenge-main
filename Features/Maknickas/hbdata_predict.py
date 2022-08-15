from __future__ import absolute_import
#from ..utils.data_utils import get_file
from os import listdir
from os.path import isfile, join
from keras.utils import np_utils
import numpy as np
import csv
import scipy.io.wavfile as wav
from numpy.lib import stride_tricks
import scipy.signal
from numpy import inf
import random
import math


def get_file_list(validation=True, test_split=0.2, seed=113):
	path = "/home/keras/examples/hbdata/"
        labels_file = path + "labels.txt"
        normal_postfix = "normal/Atraining_normal/"
        murmur_postfix = "murmur/Atraining_murmur/"

        normal_path = path+normal_postfix
        murmur_path = path+murmur_postfix

        normal_files = [f for f in listdir(normal_path) if isfile(join(normal_path, f))]
        murmur_files = [f for f in listdir(murmur_path) if isfile(join(murmur_path, f))]

        random.seed(seed)
        random.shuffle(normal_files)
        random.seed(seed)
        random.shuffle(murmur_files)

        normal_files_test = normal_files[:int(len(normal_files) * (1 - test_split))]
        normal_files_validation = normal_files[int(len(normal_files) * (1 - test_split)):]

        murmur_files_test =  murmur_files[:int(len(murmur_files) * (1 - test_split))]
        murmur_files_validation = murmur_files[int(len(murmur_files) * (1 - test_split)):]

        normal_len_test = len(normal_files_test)
        murmur_len_test = len(murmur_files_test)

        normal_len_validation = len(normal_files_validation)
        murmur_len_validation = len(murmur_files_validation)


        print("Normal files test: "+str(normal_len_test))
        print("Murmur files test: "+str(murmur_len_test))
        print("Normal files validation: "+str(normal_len_validation))
        print("Murmur files validation: "+str(murmur_len_validation))

	if validation:
		normal = normal_files_validation
		murmur = murmur_files_validation
	else:
		normal = normal_files_test
		murmur = murmur_files_test

	return (normal, murmur), (normal_path, murmur_path)	

def data_cols_from_file(filename, label, test_split=0.2, seed=113, max_frames=2, width=129, height=129):
	first_mels = mels(audiopath=filename, length=width, binsize=height)

        w, h = first_mels[0].shape

        if len(first_mels) < max_frames:
                X = np.zeros((len(first_mels), 1, w, h*3))
        else:
                X = np.zeros((max_frames, 1, w, h*3))

        labels = np.zeros((len(first_mels), 1), dtype=np.int8)
        file_num = 0

        result = first_mels
	result_delta = mels(audiopath=filename, length=width, deltas=1, binsize=height)
	result_deltadelta = mels(audiopath=filename, length=width, deltas=2, binsize=height)
#       if len(result) == 1.0:
#               X[0][0] = result[0]
#       else:
#               X[0][0] = result[1]
	for add, add_delta, add_deltadelta in zip(result, result_delta, result_deltadelta):
        	col = 0
                for q in range(w):
	                X[file_num:,:,:,col] = add[:,q]
                        X[file_num:,:,:,col+1] = add_delta[:,q]
                        X[file_num:,:,:,col+2] = add_deltadelta[:,q]
                        col += 3

                file_num += 1
                if file_num >= max_frames:
                        break
	return (X, labels)


def data_from_file(filename, label, test_split=0.2, seed=113, max_frames=2, width=129, height=129):

	first_mels = mels(audiopath=filename, length=width, binsize=height)

	w, h = first_mels[0].shape

	if len(first_mels) < max_frames:
		X = np.zeros((len(first_mels), 3, w, h))
	else:
		X = np.zeros((max_frames, 3, w, h))

	labels = np.zeros((len(first_mels), 1), dtype=np.int8)
	file_num = 0

	result = first_mels
#	if len(result) == 1.0:
#		X[0][0] = result[0]
#	else:
#		X[0][0] = result[1]
	for add in result:
		X[file_num][0] = add
		labels[file_num] = label
		file_num += 1
		if file_num >= max_frames:
			break
			

	file_num = 0
	result_delta = mels(audiopath=filename, length=width, deltas=1, binsize=height)
	
#	if len(result_delta) == 1.0:
#                X[0][1] = result[0]
#        else:
#                X[0][1] = result[1]

	for add in result_delta:
		X[file_num][1] = add
		file_num += 1
		if file_num >= max_frames:
                        break

	file_num = 0

	result_deltadelta = mels(audiopath=filename, length=width, deltas=2, binsize=height)
	
#	if len(result_deltadelta) == 1.0:
#                X[0][2] = result[0]
#        else:
#                X[0][2] = result[1]

	for add in result_deltadelta:
		X[file_num][2] = add
		file_num += 1
		if file_num >= max_frames:
                        break
              

	return (X, labels)

""" short time fourier transform of audio signal """
def stft(sig, frameSize, overlapFac=0.5, window=np.hanning):
    win = window(frameSize)
    hopSize = int(frameSize - np.floor(overlapFac * frameSize))
    
    # zeros at beginning (thus center of 1st window should be for sample nr. 0)
    samples = np.append(np.zeros(np.floor(frameSize/2.0)), sig)    
    # cols for windowing
    cols = np.ceil( (len(samples) - frameSize) / float(hopSize)) + 1
    # zeros at end (thus samples can be fully covered by frames)
    samples = np.append(samples, np.zeros(frameSize))
    
    frames = stride_tricks.as_strided(samples, shape=(cols, frameSize), strides=(samples.strides[0]*hopSize, samples.strides[0])).copy()
    frames *= win
    
    return np.fft.rfft(frames)    

def normalized(a): 
    #a -= np.mean(a, axis = 0)
    #cov = np.dot(a.T, a) / a.shape[0]
    #U,S,V = np.linalg.svd(cov)
    #Xrot = np.dot(a, U)
    #Xwhite = Xrot / np.sqrt(S + 1e-5)
    a = (a - a.min())/(a.max()-a.min())
    return a
 
""" scale frequency axis logarithmically """    
def logscale_spec(spec, sr=44100, factor=20.):
    timebins, freqbins = np.shape(spec)

    scale = np.linspace(0, 1, freqbins) ** factor
    scale *= (freqbins-1)/max(scale)
    scale = np.unique(np.round(scale))
    
    # create spectrogram with new freq bins
    newspec = np.complex128(np.zeros([timebins, len(scale)]))
    for i in range(0, len(scale)):
        if i == len(scale)-1:
            newspec[:,i] = np.sum(spec[:,scale[i]:], axis=1)
        else:        
            newspec[:,i] = np.sum(spec[:,scale[i]:scale[i+1]], axis=1)
    
    # list center freq of bins
    allfreqs = np.abs(np.fft.fftfreq(freqbins*2, 1./sr)[:freqbins+1])
    freqs = []
    for i in range(0, len(scale)):
        if i == len(scale)-1:
            freqs += [np.mean(allfreqs[scale[i]:])]
        else:
            freqs += [np.mean(allfreqs[scale[i]:scale[i+1]])]
    
    return newspec, freqs

def get_numresults(audiopath, length=1024, binsize=256):
    samplerate, samples = wav.read(audiopath)
    s = stft(samples, binsize)

    sshow, freq = logscale_spec(s, factor=1.0, sr=samplerate)

#    sshow = sshow[np.all(sshow != 0, axis=1)]

    ims = 20.*np.log10(np.abs(sshow)/10e-6) # amplitude to decibel
    ims[ims == -inf] = 0
    ims[ims == inf] = 0

#    ims = ims[np.all(ims != 0, axis=1)]

    ims = np.transpose(ims)

    freqbins, timebins = np.shape(ims)
    num_results = math.ceil(float(timebins)/float(length))

    return num_results

def mels(audiopath, length=1024, binsize=256, deltas=0):
    samplerate, samples = wav.read(audiopath)
    s = stft(samples, binsize)

    sshow, freq = logscale_spec(s, factor=1.0, sr=samplerate)

#    sshow = sshow[np.all(sshow != 0, axis=1)]

    ims = 20.*np.log10(np.abs(sshow)/10e-6) # amplitude to decibel
    ims[ims == -inf] = 0
    ims[ims == inf] = 0

#    ims = ims[np.all(ims != 0, axis=1)]

    ims = np.transpose(ims)

    if deltas > 0:
        ims = delta(ims, order=deltas)

    freqbins, timebins = np.shape(ims)

    #melss = np.zeros((freqbins, timebins)
    melss = np.array(ims)

    num_results = math.ceil(float(timebins)/float(length))

    result = np.zeros((num_results, freqbins, length))
    
    start = 0
    end = length
    for z in range(int(num_results)):
        x = 0
        for q in melss:
            if length > len(q[start:end]):
                zeros = np.zeros((length-len(q[start:end])))
    #           zeros = q[0:length-len(q)]
                q = np.concatenate((q, zeros))
            result[z][x] = q[start:end]
            x += 1
        zero_bool = result[z] == 0
        result[z][zero_bool] = 1e-5
#        result[z] = normalized(result[z])
        start += length
        end += length

#    zero_bool = result == 0
#    result[zero_bool] = 1e-5

    return normalized(result)
def delta(data, width=9, order=1, axis=-1, trim=True):
    data = np.atleast_1d(data)

    if width < 3 or np.mod(width, 2) != 1:
        print('width must be an odd integer >= 3')

    if order <= 0 or not isinstance(order, int):
        print('order must be a positive integer')

    half_length = 1 + int(width // 2)
    window = np.arange(half_length - 1., -half_length, -1.)

    # Normalize the window so we're scale-invariant
    window /= np.sum(np.abs(window)**2)

    # Pad out the data by repeating the border values (delta=0)
    padding = [(0, 0)] * data.ndim
    width = int(width)
    padding[axis] = (width, width)
    delta_x = np.pad(data, padding, mode='edge')

    for _ in range(order):
        delta_x = scipy.signal.lfilter(window, 1, delta_x, axis=axis)

    # Cut back to the original shape of the input data
    if trim:
        idx = [slice(None)] * delta_x.ndim
        idx[axis] = slice(- half_length - data.shape[axis], - half_length)
        delta_x = delta_x[idx]

    return delta_x
