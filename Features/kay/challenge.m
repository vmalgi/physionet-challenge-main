function classifyResult = challenge(recordName)
%
% Sample entry for the 2016 PhysioNet/CinC Challenge.
%
% INPUTS:
% recordName: string specifying the record name to process
%
% OUTPUTS:
% classifyResult: integer value where
%                     1 = abnormal recording
%                    -1 = normal recording
%                     0 = unsure (too noisy)
%
% To run your entry on the entire training set in a format that is
% compatible with PhysioNet's scoring enviroment, run the script
% generateValidationSet.m
%
% The challenge function requires that you have downloaded the challenge
% data 'training_set' in a subdirectory of the current directory.
%    http://physionet.org/physiobank/database/challenge/2016/
%
% This dataset is used by the generateValidationSet.m script to create
% the annotations on your training set that will be used to verify that
% your entry works properly in the PhysioNet testing environment.
%
%
% Version 1.0
%
%
% Written by: Chengyu Liu, Fubruary 21 2016
%             chengyu.liu@emory.edu
%
% Last modified by:
%
%

%% Load the trained parameter matrices for Springer's HSMM model.
% The parameters were trained using 409 heart sounds from MIT heart
% sound database, i.e., recordings a0001-a0409.
load('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');

%% Load data and resample data
springer_options   = default_Springer_HSMM_options;
[PCG, Fs1] = audioread([recordName '.wav']);  % load data
PCGResampled      = resample(PCG,springer_options.audio_Fs,Fs1); % resample to springer_options.audio_Fs (1000 Hz)

%% Running runSpringerSegmentationAlgorithm.m to obtain the assigned_states
[assignedStates] = runSpringerSegmentationAlgorithm(PCGResampled, springer_options.audio_Fs, Springer_B_matrix, Springer_pi_vector, Springer_total_obs_distribution, false); % obtain the locations for S1, systole, s2 and diastole

%Set up and get features
featuresObject = Features;
featuresObject = initialize_features(featuresObject,assignedStates,PCGResampled,springer_options.audio_Fs);
featuresObject = pre_process_PCG(featuresObject,490,2,20,2);
featuresObject = get_assignedStates_in_matrix_form(featuresObject);
featuresObject = cwt_morlet_features(featuresObject,2.4,3,7,3,7,[30,450]);
featuresObject = mfcc_features(featuresObject,21,20,40,0.1,3,7,3,7);
featuresObject = interbeat_features(featuresObject);
% featuresObject = time_domain_complexity_features(featuresObject,[2,3,4,8,20],[1]);
featuresObject = spectral_complexity_features(featuresObject,5,0.1,3,7,3,7);

features = featuresObject.featuresVector;

load('FeaturesSelectedAllTest3855.mat')
%Normalize the features
features = (features - meanOfSignal)./(standardDeviation);

%t Test to get only significant features
features = features(significant);

%Remove Covariance greater than 90% between features
features(removedFeatures) = [];

%PCA
features = (features*PCAcoeffs);


%% Running classifyFromHsIntervals.m to obtain the final classification result for the current recording
classifyResult = classifyFromHsIntervals(features);
