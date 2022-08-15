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
%

%% Load the trained parameter matrices for Springer's HSMM model.
% The parameters were trained using 409 heart sounds from MIT heart
% sound database, i.e., recordings a0001-a0409.
load('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');
load('Params.mat');
load('Classifier');
%% Load data and resample data
springer_options   = default_Springer_HSMM_options;
[ad,fs]=audioread([recordName '.wav']);
ad_resampled = resample(ad,springer_options.audio_Fs,fs); % resample to springer_options.audio_Fs (1000 Hz)
aux=runSpringerSegmentationAlgorithm(ad_resampled,springer_options.audio_Fs,Springer_B_matrix,Springer_pi_vector,Springer_total_obs_distribution,false);
ad_resampled=schmidt_spike_removal(ad_resampled,springer_options.audio_Fs);
ad_resampled=normalise_signal(ad_resampled);
F=extractFeaturesFromPCG(aux,ad_resampled);
%
% Pad with zeros to complete 172 cycles, if shorter.
st=size(F);
if st(2)<172
    Fp=F;
    Fp(:,st(2)+1:172,:)=0;   %padarray(F,[0 172-st(2) 0],0,'post');
else
    Fp=F(:,1:172,:); % cut tail if longer
end
%
%
Ft=cat(3,Fp,Fp); % duplicate to avoid recording for singleton 4-way tensor
% Reshape into 4-way tensor
for i=1:4
   T4(1:142,i,:,:)=Ft(142*(i-1)+1:142*i,:,:);
end
% Fill trailing zeros with copies of leading non-zeros
T4=fillTensor4(T4);
%
% Run classifier
%
thr=0.225;
alpha=0;
f=grabFeatures4(T4,U,indx);
f=f(1,:);
[~,prob0]=predict(Classifier,[f f.^2]);
prob = prob0(1,2);
if prob<=thr*(1-alpha)
    classifyResult=-1;
%elseif prob<=thr1*(1+alpha) & prob>thr1*(1-alpha)
%    classifyResult=0;
else
    classifyResult=1;
end
end
%

