function SQI=get_detailed_segmentation_features(audio_data,Fs)
SQI=table; 

% Variance
SQI.audio_var = get_variance_score(audio_data);

% Linear predictive coefficients
[SQI.audio_lpc10,~] = lpc(double(audio_data),10);
[SQI.audio_lpc6,  SQI.audio_lsf6] = lpc_lsf_coeff_modified(audio_data);

% Entropy based features
%[SQI.audio_shannon,SQI.audio_tsallis,SQI.audio_renyi]=get_prob_entropy(audio_data);

%% Get Power Spectral Density
p = nextpow2(Fs*2);
window_size = 2^p;
window_size=min(length(audio_data),window_size); 
[Pxx,F] = pwelch(audio_data,window_size,window_size/2,window_size,Fs);



%% Power spectral density derived features
% % Peaks features
% pxx_smooth=get_maf(Pxx);
% nb_higherPks_MAF=2;
% [SQI.power_maf_num_peaks,  SQI.power_maf_f_peaks, SQI.power_maf_diff_peaks]...
%     = peaks_features_modified(pxx_smooth,F, nb_higherPks_MAF);
% 
% [fi_tot,SQI.power_gmm_parameters]=get_gmm(Pxx,F);
% nb_higherPks_GMM=2;
% [SQI.power_gmm_num_peaks,  SQI.power_gmm_f_peaks, SQI.power_gmm_diff_peaks]...
%     = peaks_features_modified(fi_tot,F, nb_higherPks_GMM);

% Linear regression line
%[SQI.power_regression_slope,SQI.power_regression_intercept,SQI.power_regression_r2]=get_regression_line(Pxx,F); 

% Spectral composition
[SQI.power_b1, SQI.power_b2, SQI.power_b3, SQI.power_b4, SQI.power_b5, ...
    SQI.power_b6, SQI.power_b7, SQI.power_b8, SQI.power_b9, SQI.power_b10] = get_spectral_composition(Pxx,F);
[SQI.power_f_mean,SQI.power_std,SQI.power_f_med,SQI.power_bw,SQI.power_f_p25,...
    SQI.power_f_p75,SQI.power_f_IQR,SQI.power_000_100,SQI.power_100_200,...
    SQI.power_200_400,SQI.power_400_600,SQI.power_600_800,SQI.power_800_1000,...
    SQI.power_1000_1200,SQI.power_1200_1400,SQI.power_1400_1600,SQI.power_1600_1800,...
    SQI.power_1800_2000,SQI.power_tp]=get_spectrum_parameters(Pxx,F);   

% % Modified power spectral density centroid
% SQI.power_freq_centroid=get_power_freq_centroid(Pxx,F); 
% 
% % Spectral entropy
% SQI.power_hn = get_spectral_entropy(Pxx);

% Dominant frequency freatures
[SQI.power_max_pow, SQI.power_max_freq,SQI.power_ratio_max] = get_dominant_frequency_features(Pxx,F);

% Variance
SQI.power_var = get_variance_score(Pxx);
end