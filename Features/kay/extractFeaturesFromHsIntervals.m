function features = extractFeaturesFromHsIntervals(assignedStates, PCG)


% This function calculate features based on the segmented heart signal
% and the PCG
%
% INPUTS:
% assignedStates: the array of state values assigned to the sound recording.
% PCG: resampled sound recording with 1000 Hz.
%
% OUTPUTS:
% features: the obtained 20 features for the current sound recording
%
%
% Written by: Chengyu Liu, January 22 2016
%             chengyu.liu@emory.edu
%
% Last modified by: Edmund Kay
%                   ek360@cam.ac.uk



%% %%%%%%%%%%%%%%%%%% DEFAULT inputs for sample frequencies %%%%%%%%%%%%%%%%%



springer_options = default_Springer_HSMM_options;
Fs = springer_options.audio_Fs;
downsampleFs = springer_options.audio_segmentation_Fs;
%[resampled_assignedStates] = integer_resample(assignedStates,...
% springer_options.audio_segmentation_Fs,springer_options.audio_Fs);



%% %%%%%%%%%%%%%%%%%% Get assigned states into matrix form %%%%%%%%%%%%%%%%%%



indx = find(abs(diff(assignedStates))>0); % find the locations with changed states
switch assignedStates(indx(1)+1)
    case 4
        K=1;
    case 3
        K=2;
    case 2
        K=3;
    case 1
        K=0;
end
K=K+1;

indx2                = indx(K:end);
rem                  = mod(length(indx2),4);
indx2(end-rem+1:end) = [];
A                    = reshape(indx2,4,length(indx2)/4)';
A = A+1;



%% %%%%%%%%%%%%%%%%%% Filter PCG and remove spikes %%%%%%%%%%%%%%%%%%%%%%%%%%%



PCG = butterworth_low_pass_filter(PCG,2,490,Fs, false);
PCG = butterworth_high_pass_filter(PCG,2,20,Fs);
PCG = schmidt_spike_removal(PCG,Fs);



%% %%%%%%%%%%%%%%%%%% Pure Wavelet transform for representation %%%%%%%%%%%%%%

wavelet_name ='morl';
freqrange = [30 450];
voicesPerOctave = 3;
scales = scalesFromFrequencies(wavelet_name,freqrange,voicesPerOctave,Fs);
% Audio needs to be longer than 1 second
if(length(PCG)< Fs*1.025)
    PCG = [PCG; zeros(round(0.025*Fs),1)];
end
[cd] = cwt(PCG,scales,wavelet_name);
cd = cd';
TFR = abs(cd);
normalizedTFR = TFR;
normalizedTFR = per_frequency_normalization(TFR);

NS1 = 3;
NSystole = 7;
NS2 = 3;
NDiastole = 7;

[S1Segments] = get_timings_of_selected_segments(A(:,1:2),NS1);
[SystoleSegments] = get_timings_of_selected_segments(A(:,2:3),NSystole);
[S2Segments] = get_timings_of_selected_segments(A(:,3:4),NS2);
[DiastoleSegments] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastole);

meanCWTCell = cell(size(normalizedTFR,2),1);
meanCWT = [];
for i = 1:size(normalizedTFR,2)
    meanCWTCell{i} = get_feature_over_beat(normalizedTFR(:,i),S1Segments,...
                                SystoleSegments,S2Segments,DiastoleSegments);
    meanCWT = [meanCWT,meanCWTCell{i}];
end



%% %%%%%%%%%%%%%%%%%%%% Estimation of Spectral Flux from CWT %%%%%%%%%%%%%%%%%%%%

% wavelet_name ='morl';
% freqrange = [30 450];
% voicesPerOctave = 3;
% scales = scalesFromFrequencies(wavelet_name,freqrange,voicesPerOctave,Fs);
% % Audio needs to be longer than 1 second
% if(length(PCG)< Fs*1.025)
%     PCG = [PCG; zeros(round(0.025*Fs),1)];
% end
% [cd] = cwt(PCG,scales,wavelet_name);
% cd = cd';
% TFR = abs(cd);
% normalizedTFR = TFR;
% normalizedTFR = per_frequency_normalization(TFR);
% 
% NS1 = 10;
% NSystole = 30;
% NS2 = 30;
% NDiastole = 30;
% 
% [S1Segments] = get_timings_of_selected_segments(A(:,1:2),NS1);
% [SystoleSegments] = get_timings_of_selected_segments(A(:,2:3),NSystole);
% [S2Segments] = get_timings_of_selected_segments(A(:,3:4),NS2);
% [DiastoleSegments] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastole);
% 
% meanCWTCell = cell(size(normalizedTFR,2),1);
% meanCWTForFlux = [];
% for i = 1:size(normalizedTFR,2)
%     meanCWTCell{i} = get_feature_over_beat(normalizedTFR(:,i),S1Segments,...
%                                 SystoleSegments,S2Segments,DiastoleSegments);
%     meanCWTForFlux = [meanCWTForFlux,meanCWTCell{i}];
% end
% figure
% plot(meanCWTForFlux)
% 'a'


%% %%%%%%%%%%%%%%%%%% Homomorphic Envelope for deriving features %%%%%%%%%%%%%

% 
% homomorphicEnvelope = Homomorphic_Envelope_with_Hilbert(PCG, Fs);
% homomorphicEnvelope = normalise_signal(homomorphicEnvelope);
% 
% NS1Envelope = 5;
% NSystoleEnvelope = 5;
% NS2Envelope = 100;
% NDiastoleEnvelope = 5;
% 
% [S1Segments] = get_timings_of_selected_segments(A(:,1:2),NS1Envelope);
% [SystoleSegments] = get_timings_of_selected_segments(A(:,2:3),NSystoleEnvelope);
% [S2Segments] = get_timings_of_selected_segments(A(:,3:4),NS2Envelope);
% [DiastoleSegments] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastoleEnvelope);
% 
% meanEnvelope = get_feature_over_beat(homomorphicEnvelope,...
%     S1Segments,SystoleSegments,S2Segments,DiastoleSegments);
% 
% meanEnvelope = meanEnvelope./max(meanEnvelope);

  
% meanFHSAmplitude = mean([mean(meanEnvelope(1:NS1Envelope)),...
%    mean(meanEnvelope(1+NS1Envelope+NSystoleEnvelope:NS1Envelope+NSystoleEnvelope+NS2Envelope))]);
% % meanS1Amplitude = mean(meanEnvelope(1:NS1Envelope));
% medSystolicAmplitude = meanEnvelope(ceil(NS1Envelope+(NSystoleEnvelope/2)));
% % meanS2Amplitude = mean(meanEnvelope(1+NS1Envelope+NSystoleEnvelope:NS1Envelope+NSystoleEnvelope+NS2Envelope));
% medDiastolicAmplitude = meanEnvelope(ceil(NS1Envelope+NSystoleEnvelope+NS2Envelope+(NDiastoleEnvelope/2)));


% sdSystolicAmplitude = std(meanEnvelope(NS1Envelope+1:NS1Envelope+NSystoleEnvelope));
% sdDiastolicAmplitude = std(meanEnvelope(1+NS1Envelope+NSystoleEnvelope+NS2Envelope:end));


% figure
% plot(meanEnvelope)

% envelopeFeatures = [sdSystolicAmplitude,sdDiastolicAmplitude];

%% %%%%%%%% High Resolution Wavelet Transform for deriving features %%%%%%%%%%

% 
% 
% wavelet_name ='morl';
% freqrange = [30 450];
% voicesPerOctave = 12;
% scales = scalesFromFrequencies(wavelet_name,freqrange,voicesPerOctave,Fs);
% freqs = scal2frq(scales,wavelet_name,1./Fs);
% % Audio needs to be longer than 1 second
% if(length(PCG)< Fs*1.025)
%     PCG = [PCG; zeros(round(0.025*Fs),1)];
% end
% [cd] = cwt(PCG,scales,wavelet_name);
% cd = cd';
% TFR = abs(cd);
% 
% normalizedTFR = TFR;
% % normalizedTFR = per_frequency_normalization(TFR);
% 
% NS1HD = 5;
% NSystoleHD = 15;
% NS2HD = 5;
% NDiastoleHD = 15;
% 
% [S1Segments] = get_timings_of_selected_segments(A(:,1:2),NS1HD);
% [SystoleSegments] = get_timings_of_selected_segments(A(:,2:3),NSystoleHD);
% [S2Segments] = get_timings_of_selected_segments(A(:,3:4),NS2HD);
% [DiastoleSegments] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastoleHD);
% 
% meanCWTCellHD = cell(size(normalizedTFR,2),1);
% for i = 1:size(normalizedTFR,2)
%     meanCWTCellHD{i} = get_feature_over_beat(normalizedTFR(:,i),S1Segments,...
%                                 SystoleSegments,S2Segments,DiastoleSegments);
% end



%% %%%%%%%%% Derive Features from HD wavelet transform %%%%%%%%%%%%%%%%%%%%%%%%


% 
% S1Features = zeros(NS1HD,length(meanCWTCellHD));
% SystolicFeatures = zeros(NSystoleHD,length(meanCWTCellHD));
% S2Features = zeros(NS2HD,length(meanCWTCellHD));
% DiastolicFeatures = zeros(NDiastoleHD,length(meanCWTCellHD));
% 
% for i = 1:length(meanCWTCellHD)
%     S1Features(:,i) = meanCWTCellHD{i}(1:NS1HD);
%     SystolicFeatures(:,i) = meanCWTCellHD{i}(NS1HD+1:NS1HD+NSystoleHD);
%     S2Features(:,i) = meanCWTCellHD{i}(NS1HD+NSystoleHD+1:NS1HD+NSystoleHD+NS2HD);
%     DiastolicFeatures(:,i) = meanCWTCellHD{i}(NS1HD+NSystoleHD+NS2HD+1:NS1HD+NSystoleHD+NS2HD+NDiastoleHD);
% end
% 
% % S1Freq = fundamental_frequency(S1Features,freqs);
% % S1Bandwidth = bandwidth_frequency(S1Features,freqs,2,NS1HD);
% 
% S2Freq = fundamental_frequency(S2Features,freqs);
% % S2Bandwidth = bandwidth_frequency(S2Features,freqs,2,NS2HD);
% 
% SystoleFreq = fundamental_frequency(SystolicFeatures,freqs);
% % SystoleBandwidth = bandwidth_frequency(SystolicFeatures,freqs,2,NSystoleHD);
% positionMaxFreqSystole = position_max_freq(SystolicFeatures,freqs,NSystoleHD);
% 
% DiastoleFreq = fundamental_frequency(DiastolicFeatures,freqs);
% % DiastoleBandwidth = bandwidth_frequency(DiastolicFeatures,freqs,2,NDiastoleHD);
% positionMaxFreqDiastole = position_max_freq(DiastolicFeatures,freqs,NDiastoleHD);
% 
% SysS2Freq = SystoleFreq/S2Freq;
% DiaS2Freq = DiastoleFreq/S2Freq;
% 
% derivedCWTFeatures = [positionMaxFreqSystole,positionMaxFreqDiastole,...
%                        SysS2Freq,DiaS2Freq];
% 
            
                   
                   
%% %%%%%%%%%%%%%%%%%%%%%%%% Interbeat Features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m_RR        = round(mean(diff(A(:,1))));             % mean value of RR intervals
sd_RR       = round(std(diff(A(:,1))));              % standard deviation (SD) value of RR intervals
mean_IntS1  = round(mean(A(:,2)-A(:,1)));            % mean value of S1 intervals
sd_IntS1    = round(std(A(:,2)-A(:,1)));             % SD value of S1 intervals
%mean_IntS2  = round(mean(A(:,4)-A(:,3)));            % mean value of S2 intervals
sd_IntS2    = round(std(A(:,4)-A(:,3)));             % SD value of S2 intervals
mean_IntSys = round(mean(A(:,3)-A(:,2)));            % mean value of systole intervals
sd_IntSys   = round(std(A(:,3)-A(:,2)));             % SD value of systole intervals
mean_IntDia = round(mean(A(2:end,1)-A(1:end-1,4)));  % mean value of diastole intervals
sd_IntDia   = round(std(A(2:end,1)-A(1:end-1,4)));   % SD value of diastole intervals

for i=1:size(A,1)-1
    R_SysRR(i)  = (A(i,3)-A(i,2))/(A(i+1,1)-A(i,1));
    R_DiaRR(i)  = (A(i+1,1)-A(i,4))/(A(i+1,1)-A(i,1));
    R_SysDia(i) = R_SysRR(i)/R_DiaRR(i);
    
    P_S1(i)     = sum(abs(PCG(A(i,1):A(i,2))))/(A(i,2)-A(i,1));
    P_Sys(i)    = sum(abs(PCG(A(i,2):A(i,3))))/(A(i,3)-A(i,2));
    P_S2(i)     = sum(abs(PCG(A(i,3):A(i,4))))/(A(i,4)-A(i,3));
    P_Dia(i)    = sum(abs(PCG(A(i,4):A(i+1,1))))/(A(i+1,1)-A(i,4));
    if P_S1(i)>0
        P_SysS1(i) = P_Sys(i)/P_S1(i);
    else
        P_SysS1(i) = 0;
    end
    if P_S2(i)>0
        P_DiaS2(i) = P_Dia(i)/P_S2(i);
    else
        P_DiaS2(i) = 0;
    end
end

% m_Ratio_SysRR   = mean(R_SysRR);  
sd_Ratio_SysRR  = std(R_SysRR);  
% m_Ratio_DiaRR   = mean(R_DiaRR);  
sd_Ratio_DiaRR  = std(R_DiaRR);   
% m_Ratio_SysDia  = mean(R_SysDia);
sd_Ratio_SysDia = std(R_SysDia); 

indx_sys = find(P_SysS1>0 & P_SysS1<100);   % avoid the flat line signal
if length(indx_sys)>1
    m_Amp_SysS1  = mean(P_SysS1(indx_sys)); 
    sd_Amp_SysS1 = std(P_SysS1(indx_sys));  
else
    m_Amp_SysS1  = 0;
    sd_Amp_SysS1 = 0;
end
indx_dia = find(P_DiaS2>0 & P_DiaS2<100);
if length(indx_dia)>1
    m_Amp_DiaS2  = mean(P_DiaS2(indx_dia)); 
    sd_Amp_DiaS2 = std(P_DiaS2(indx_dia)); 
else
    m_Amp_DiaS2  = 0;
    sd_Amp_DiaS2 = 0;
end

interBeatFeatures = [m_RR, sd_RR, mean_IntS1, sd_IntS1, sd_IntS2,...
    mean_IntSys, sd_IntSys, mean_IntDia, sd_IntDia, sd_Ratio_SysRR,...
    sd_Ratio_DiaRR, sd_Ratio_SysDia,m_Amp_SysS1,sd_Amp_SysS1,...
    m_Amp_DiaS2, sd_Amp_DiaS2];






%% %%%%%%%%%%%%%%% Different Time Frequency Representation %%%%%%%%%%%%%%%%%%



% t = [1:length(PCG)];
% fBins = 8;
% tempFBins = ceil(fBins./((2.*(freqrange(2)-freqrange(1)))./(Fs)));
% [TFR t freqs] = tfrpwv(PCG,t,tempFBins);
% freqs = freqs.*Fs;
% index = find(freqs>=freqrange(1) & freqs<=freqrange(2));
% TFR = TFR([index(1)-1:index(end)+1],:);
% TFR = abs(TFR)';
% freqs = freqs([index(1)-1:index(end)+1]);



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Get MFCC %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nCepstralFeatures = 21;
nFilterBanks = 20;
nFrequencyBins = 40;
fractionOverlap = 0.1;


NS1Cep = 3;
NSystoleCep = 7;
NS2Cep = 3;
NDiastoleCep = 7;

[S1SegmentsCep] = get_timings_of_selected_segments(A(:,1:2),NS1Cep);
[SystoleSegmentsCep] = get_timings_of_selected_segments(A(:,2:3),NSystoleCep);
[S2SegmentsCep] = get_timings_of_selected_segments(A(:,3:4),NS2Cep);
[DiastoleSegmentsCep] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastoleCep);

[Pxx,f] = powerSpectrumMFCC(PCG,S1SegmentsCep,SystoleSegmentsCep,S2SegmentsCep,DiastoleSegmentsCep,fractionOverlap,nFrequencyBins,Fs);
% figure
% plot(reshape(log(Pxx'+eps),[size(Pxx,1)*size(Pxx,2),1]))

H = melSpacedFilterBank(Pxx,f,nFilterBanks);
PxxFiltered = (H*Pxx)';

Cepstrum = zeros(size(PxxFiltered,1),nCepstralFeatures);
for i = 1:nCepstralFeatures
Cepstrum(:,i) = log(PxxFiltered+eps)*cos(i.*([1:size(PxxFiltered,2)]' - 0.5).*(pi/nCepstralFeatures));
end
Cepstrum = Cepstrum(:,1:end-1);

% figure 
% surf(Cepstrum','edgecolor','none')
% view([0,90])


%In order to calculate deltas we need two values before and two values
%after the cepstrum assume that heartbeat circular so append first two
%values and prepend last two values
% 
% N = 2;
% 
% CepstrumForDeltas = [Cepstrum(end-(N-1):end,:);
%     Cepstrum;
%     Cepstrum(1:N,:)];
% 
% CepstrumForDeltaDeltas = [Cepstrum(end-(N^2-1):end,:);
%     Cepstrum;
%     Cepstrum(1:N^2,:)];
% 
%                  
% deltas = zeros(size(Cepstrum));          
% for jSum = 1:N
%    deltas = deltas + ((jSum.*(CepstrumForDeltas(1+N+jSum:end-N+jSum,:) - CepstrumForDeltas(1+N-jSum:end-N-jSum,:)))./(2.*(sum([1:N].^2))));
% end
% 
% deltasForDeltaDeltas = zeros(size(CepstrumForDeltas));   
% for jSum = 1:N
%    deltasForDeltaDeltas = deltasForDeltaDeltas + ((jSum.*(CepstrumForDeltaDeltas(1+N+jSum:end-N+jSum,:) - CepstrumForDeltaDeltas(1+N-jSum:end-N-jSum,:)))./(2.*(sum([1:N].^2))));
% end
%    
% deltaDeltas = zeros(size(Cepstrum));              
% for jSum = 1:N
%    deltaDeltas = deltaDeltas + ((jSum.*(deltasForDeltaDeltas(1+N+jSum:end-N+jSum,:) - deltasForDeltaDeltas(1+N-jSum:end-N-jSum,:)))./(2.*(sum([1:N].^2))));
% end
%                  


Cepstrum = reshape(Cepstrum,[size(Cepstrum,1)*size(Cepstrum,2),1])';
% deltas = reshape(deltas,[size(deltas,1)*size(deltas,2),1])';
% deltaDeltas = reshape(deltaDeltas,[size(deltaDeltas,1)*size(deltaDeltas,2),1])';





%% %%%%%%%%%%%%%%%%%%%%%%%%%% Signal Complexity features %%%%%%%%%%%%%%%%%%%%%

% These signal complextiy features come from Schmidts 2015 paper "Acoustic
% Detection of coronary artery disease"





%% Get X vector 
m = [2,3,4,8,20];
tau = [1,2,4,6,8,10,12];


middleBeat = round(size(A,1)/2);
PCGMiddleBeat = PCG(A(middleBeat,1):A(middleBeat+1,1));
PCGForComplexity = PCGMiddleBeat(1:2:end);

trajectoryMatrixCell = cell(length(m),length(tau));
simplicityMatrix = zeros(length(m),length(tau));

for i=1:length(m)
    for j = 1:length(tau)
        p = length(PCGForComplexity) - (m(i)-1)*tau(j);
        trajectoryMatrix = zeros(p,m(i));
        for k = 1:p
            for l = 1:m(i)
                trajectoryMatrix(k,l) = PCGForComplexity(k+(l-1)*tau(j)) ;
            end
        end
        trajectoryMatrix = trajectoryMatrix./(size(trajectoryMatrix,1));
        trajectoryMatrixCell{i,j} = trajectoryMatrix;
        
        
        ComplexityMatrix = transpose(trajectoryMatrix)*trajectoryMatrix;
        lambdas = eig(ComplexityMatrix);
        lambdaSort = sort(lambdas,'descend');
        lambdaSort = lambdaSort./(sum(lambdaSort));
        entropy = -sum(lambdaSort.*log(lambdaSort+eps));
        omega = 2.^entropy;
        simplicity = 1./omega;
        simplicityMatrix(i,j) = simplicity;
        
    end
end

simplicityVector = reshape(simplicityMatrix,[size( simplicityMatrix,1)*size( simplicityMatrix,2),1])';


%     Calculate sample entropy, only works for m = 2
threshold = 3*10^-5;
phi = zeros(2,length(tau));
for i = 1:length(m)
    for j = 1:length(tau)
        trajectoryMatrix = trajectoryMatrixCell{i,j};
        C = zeros(size(trajectoryMatrix,1),1);
        for k = 1:size(trajectoryMatrix,1)
            vector1 = trajectoryMatrix(k,:);
            differenceMatrix = trajectoryMatrix - repmat(vector1,[size(trajectoryMatrix,1), 1]);
            differenceMatrixMax = max(abs(differenceMatrix),[],2);
            count = sum(differenceMatrixMax <= threshold);
            C(k) = count./size(trajectoryMatrix,1);
        end
        phi(i,j) = (1./size(trajectoryMatrix,1)).*(sum(log(C+eps)));
    end
end

approximateEntropy = zeros(size(phi,1)-1,size(phi,2));
for i = 1:size(approximateEntropy,1)
   approximateEntropy(i,:) = phi(i,:) - phi(i+1,:); 
end

approximateEntropyVector = reshape(approximateEntropy,[size( approximateEntropy,1)*size( approximateEntropy,2),1])';

% % wavelet_name ='morl';
% % freqrange = [30 450];
% % voicesPerOctave = 2;
% % scales = scalesFromFrequencies(wavelet_name,freqrange,voicesPerOctave,Fs);
% % % Audio needs to be longer than 1 second
% % if(length(PCG)< Fs*1.025)
% %     PCG = [PCG; zeros(round(0.025*Fs),1)];
% % end
% % [cd] = cwt(PCG,scales,wavelet_name);
% % cd = cd';
% % TFR = abs(cd);
% % normalizedTFR = TFR;
% % normalizedTFR = bulk_normalization(TFR);
% % 
% % NS1 = 20;
% % NSystole = 80;
% % NS2 = 20;
% % NDiastole = 80;
% % 
% % [S1Segments] = get_timings_of_selected_segments(A(:,1:2),NS1);
% % [SystoleSegments] = get_timings_of_selected_segments(A(:,2:3),NSystole);
% % [S2Segments] = get_timings_of_selected_segments(A(:,3:4),NS2);
% % [DiastoleSegments] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastole);
% % 
% % meanCWTCell = cell(size(normalizedTFR,2),1);
% % meanCWTForComplexity = [];
% % for i = 1:size(normalizedTFR,2)
% %     meanCWTCell{i} = get_feature_over_beat(normalizedTFR(:,i),S1Segments,...
% %                                 SystoleSegments,S2Segments,DiastoleSegments);
% %     meanCWTForComplexity = [meanCWTForComplexity;meanCWTCell{i}];
% % end
% % 
% % meanCWTForComplexity = (meanCWTForComplexity-min(min(meanCWTForComplexity)))./(max(max(meanCWTForComplexity))-min(min(meanCWTForComplexity)));
% 
% % figure
% % plot(reshape(meanCWTForComplexity',[size(meanCWTForComplexity,1)*size(meanCWTForComplexity,2),1])')
% 
% % mobility = std(PCGMiddleBeatPrime)./std(PCGMiddleBeat);
% % complexity = std(PCGMiddleBeatPrimePrime)./std(PCGMiddleBeatPrime);
% 
% 
% % SpecS = -sum(meanCWTForComplexity.*log(meanCWTForComplexity),2)';
% % StandardDeviation = std(meanCWTForComplexity,0,2)';
% % Skewness = skewness(meanCWTForComplexity,0,2)';
% % Kurtosis = kurtosis(meanCWTForComplexity,0,2)';
% % SystolicStandardDeviation = std(meanCWTForComplexity(:,(NS1+1):(NS1+NSystole)),0,2)';
% % SystolicSkewness = skewness(meanCWTForComplexity(:,(NS1+1):(NS1+NSystole)),0,2)';
% % SystolicKurtosis = kurtosis(meanCWTForComplexity(:,(NS1+1):(NS1+NSystole)),0,2)';
% % DiastlicStandardDeviation = std(meanCWTForComplexity(:,(NS1+NSystole+NS2+1):(NS1+NSystole+NS2+NDiastole)),0,2)';
% % DiasoticSkewness = skewness(meanCWTForComplexity(:,(NS1+NSystole+NS2+1):(NS1+NSystole+NS2+NDiastole)),0,2)';
% % DiastolicKurtosis = kurtosis(meanCWTForComplexity(:,(NS1+NSystole+NS2+1):(NS1+NSystole+NS2+NDiastole)),0,2)';
% % 
% % SpectralComplexityFeatures = [SpecS,StandardDeviation,Skewness,Kurtosis,SystolicStandardDeviation,SystolicSkewness,SystolicKurtosis,DiastlicStandardDeviation,DiasoticSkewness,DiastolicKurtosis];
% 
%For features that require only time domain signal, take middle beat
middleBeat = round(size(A,1)/2);
PCGMiddleBeat = PCG(A(middleBeat,1):A(middleBeat,4));
PCGMiddleBeatPrime = diff(PCGMiddleBeat);
% PCGMiddleBeatPrimePrime = diff(PCGMiddleBeatPrime);

% zeroCrossingPoints = 0;
% for i = 1:(length(PCGMiddleBeat)-1)
%     if PCGMiddleBeat(i)*PCGMiddleBeat(i+1) < 0
%         zeroCrossingPoints = zeroCrossingPoints+1;
%     end
% end
turningPoints = 0;
for i = 1:(length(PCGMiddleBeatPrime)-1)
    if PCGMiddleBeatPrime(i)*PCGMiddleBeatPrime(i+1) < 0
        turningPoints = turningPoints+1;
    end
end




%turning points is the only statistically significant one
simpleComplexityFeatures = [turningPoints,simplicityVector,approximateEntropyVector];




%Spectral entropy

nFrequencyBins = 5;
fractionOverlap = 0.1;


NS1SpecS = 3;
NSystoleSpecS = 7;
NS2SpecS = 3;
NDiastoleSpecS = 7;

[S1SegmentsSpecS] = get_timings_of_selected_segments(A(:,1:2),NS1SpecS);
[SystoleSegmentsSpecS] = get_timings_of_selected_segments(A(:,2:3),NSystoleSpecS);
[S2SegmentsSpecS] = get_timings_of_selected_segments(A(:,3:4),NS2SpecS);
[DiastoleSegmentsSpecS] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastoleSpecS);

[PxxSpecS,f] = powerSpectrumMFCC(PCG,S1SegmentsSpecS,SystoleSegmentsSpecS,S2SegmentsSpecS,DiastoleSegmentsSpecS,fractionOverlap,nFrequencyBins,Fs);

PxxSpecS = PxxSpecS./(max(max(PxxSpecS)));

SpecS = -sum(PxxSpecS.*log(PxxSpecS),1);
StandardDeviation = std(PxxSpecS,0,2)';
Skewness = skewness(PxxSpecS,0,2)';
Kurtosis = kurtosis(PxxSpecS,0,2)';

SpectralComplexityFeatures = [SpecS,StandardDeviation,Skewness,Kurtosis];

%% %%%%%%%%%%%%%%%%%%%%%%%%%%% LPC Features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The implementation of this comes from Hermansky 1989 (not completed as I'm
% not convinced that it will add anything to the CWT and MFCC list of
% features
% 
% nFrequencyBins = 40;
% fractionOverlap = 0.1;
% 
% % 1) Get the Power Spectrum in the same way as for MFCC 
% NS1LPC = 3;
% NSystoleLPC = 7;
% NS2LPC = 3;
% NDiastoleLPC = 7;
% 
% [S1SegmentsLPC] = get_timings_of_selected_segments(A(:,1:2),NS1LPC);
% [SystoleSegmentsLPC] = get_timings_of_selected_segments(A(:,2:3),NSystoleLPC);
% [S2SegmentsLPC] = get_timings_of_selected_segments(A(:,3:4),NS2LPC);
% [DiastoleSegmentsLPC] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastoleLPC);
% 
% [PxxRaw,fRaw] = powerSpectrumMFCC(PCG,S1SegmentsLPC,SystoleSegmentsLPC,S2SegmentsLPC,DiastoleSegmentsLPC,fractionOverlap,nFrequencyBins,Fs);
% 
% % 2) Warp Pxx along frequency axis
% [Pxx,f] = bark_frequency_warp(PxxRaw,fRaw);
% % 3) Convolve warped Pxx with critical band masking curve
% [theta, f] = critical_band_convolution(Pxx,f);
% % 4) Preempasize equal loudness
% [phi] = equal_loudness_preemphasis(theta,f);
% % 5) Intensity loudness power law
% finalPxx = phi.^(0.33);


%% %%%%%%%%%%%%%%%%%%%%%% Final Feature Selection %%%%%%%%%%%%%%%%%%%%%%%%%%%%



% features = [meanCWT,envelopeFeatures,derivedCWTFeatures,interBeatFeatures];
features = [meanCWT,Cepstrum,interBeatFeatures,simpleComplexityFeatures,SpectralComplexityFeatures];
% features = [SpectralComplexityFeatures];
% features = [Cepstrum,deltas,deltaDeltas];
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plot features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% figure
% plot(features)
% 'a'
% ylims = get(gca,'ylim');
% hold on
% plot([20 20],ylims,'k--')
% hold on
% plot([40 40],ylims,'k--')
% hold on
% plot([60 60],ylims,'k--')
% hold on
% plot([80 80],ylims,'k--')
% hold on
% plot([100 100],ylims,'k--')
% hold on
% plot([120 120],ylims,'k--')
% hold on
% plot([140 140],ylims,'k--')
% hold on
% plot([160 160],ylims,'k--')
% hold on
% plot([180 180],ylims,'k--')
% hold on
% plot([200 200],ylims,'k--')
% hold on
% plot([220 220],ylims,'k--')
% hold on
% plot([227 227],ylims,'k--')
% hold on
% plot([234 234],ylims,'k--')
% hold on
% plot([241 241],ylims,'k--')
% hold on
% plot([248 248],ylims,'k--')
% hold on
% plot([255 255],ylims,'k--')
% hold on
% plot([262 262],ylims,'k--')
% hold on
% plot([269 269],ylims,'k--')
% hold on
% plot([276 276],ylims,'k--')
% hold on


