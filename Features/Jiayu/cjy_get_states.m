gload('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');
springer_options    = default_Springer_HSMM_options;
N = length(signal);
features = zeros(N,20);
for i = 1 :N
[assigned_states{i}] =  runSpringerSegmentationAlgorithm(signal{i}, springer_options.audio_Fs, Springer_B_matrix, Springer_pi_vector, Springer_total_obs_distribution, false); % obtain the locations for S1, systole, s2 and diastole
features(i,:)  = extractFeaturesFromHsIntervals(assigned_states{i},signal{i});
end
