classdef Features
    properties
        featuresVector
        assignedStates
        PCG
        Fs
    end
    methods
        
        function obj = initialize_features(obj,assignedStates,PCG,Fs)
            obj.featuresVector = [];
            obj.assignedStates = assignedStates;
            obj.PCG = PCG;
            obj.Fs = Fs;
        end
        
        function obj = pre_process_PCG(obj,LowPassFreq,FilterOrderLowPass,HighPassFreq,FilterOrderHighPass)
            obj.PCG = Features.butterworth_low_pass_filter(obj.PCG,FilterOrderLowPass,LowPassFreq,obj.Fs);
            obj.PCG = Features.butterworth_high_pass_filter(obj.PCG,FilterOrderHighPass,HighPassFreq,obj.Fs);
            obj.PCG = Features.schmidt_spike_removal(obj.PCG,obj.Fs);
        end
        
        function obj = get_assignedStates_in_matrix_form(obj)
            indx = find(abs(diff(obj.assignedStates))>0); % find the locations with changed states
            switch obj.assignedStates(indx(1)+1)
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
            
            obj.assignedStates = A;
        end
        
        
        %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Continuous Wavelet Transform %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = cwt_morlet_features(obj,voicesPerOctave,NS1,NSystole,NS2,NDiastole,freqRange)
            
            wavelet_name = 'morl';
            scales = scales_from_frequencies(wavelet_name,freqRange,voicesPerOctave,obj.Fs);
            if(length(obj.PCG)< obj.Fs*1.025)
                obj.PCG = [obj.PCG; zeros(round(0.025*obj.Fs),1)];
            end
            [cd] = cwt(obj.PCG,scales,wavelet_name);
            cd = cd';
            TFR = abs(cd);
            normalizedTFR = TFR;
            normalizedTFR = Features.per_frequency_normalization(TFR);
            
            NTot = NS1 + NSystole + NS2 + NDiastole;
            
            [S1Segments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,1:2),NS1);
            [SystoleSegments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,2:3),NSystole);
            [S2Segments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,3:4),NS2);
            [DiastoleSegments] = Features.get_timings_of_selected_segments([obj.assignedStates(1:end-1,4) obj.assignedStates(2:end,1)],NDiastole);
            
            meanCWTAllBeats = zeros(size(obj.assignedStates,1)-1,NTot,size(normalizedTFR,2));
            meanCWTMatrix = zeros(NTot,size(normalizedTFR,2));
            for i = 1:size(normalizedTFR,2)
                meanCWTAllBeats(:,:,i) = Features.get_feature_over_beat(normalizedTFR(:,i),S1Segments,...
                    SystoleSegments,S2Segments,DiastoleSegments);
                
                meanCWTMatrix(:,i) = Features.select_correlated_beats_and_mean(meanCWTAllBeats(:,:,i),1);
            end
            
            % meanCWTMatrix = select_correlated_beats_and_mean(meanCWTAllBeats,1);
            meanCWT = reshape(meanCWTMatrix,[1,size(meanCWTMatrix,1)*size(meanCWTMatrix,2)]);
            
            obj.featuresVector = [obj.featuresVector,meanCWT];
            
            function scales = scales_from_frequencies(wavelet,freqrange,voicesPerOctave,fs)
                
                cfreq = centfrq(wavelet);
                
                minscale = cfreq./(freqrange(2).*(1./fs));
                maxscale = cfreq./(freqrange(1).*(1./fs));
                
                tempMinScale = floor(voicesPerOctave*log2(minscale));
                tempMaxScale = ceil(voicesPerOctave*log2(maxscale));
                
                scales = ((2^(1/voicesPerOctave)).^(tempMinScale:tempMaxScale));
                
            end
            
            
        end
        
        %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Mel Frequency Cepstral Coeffs %%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = mfcc_features(obj,nCepstralFeatures,nFilterBanks,nFrequencyBins,fractionOverlap,NS1,NSystole,NS2,NDiastole)
            [S1Segments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,1:2),NS1);
            [SystoleSegments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,2:3),NSystole);
            [S2Segments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,3:4),NS2);
            [DiastoleSegments] = Features.get_timings_of_selected_segments([obj.assignedStates(1:end-1,4) obj.assignedStates(2:end,1)],NDiastole);
            
            [Pxx,f] = Features.power_spectrum_MFCC(obj.PCG,S1Segments,SystoleSegments,S2Segments,DiastoleSegments,fractionOverlap,nFrequencyBins,obj.Fs);
            
            H = mel_spaced_filter_bank(f,nFilterBanks);
            
%             figure
%             for i = 1:size(H,1)
%                 plot(f,H(i,:))
%                 xlabel('Frequency/ Hz')
%                 set(gca,'XTick',[0,100,200,300,400,500])
%                 set(gca,'YTick',[0,0.2,0.4,0.6,0.8,1])
%                 hold on
%             end
            
            
            PxxFiltered = (H*Pxx)';
            
            Cepstrum = zeros(size(PxxFiltered,1),nCepstralFeatures);
            for jj = 1:nCepstralFeatures
                Cepstrum(:,jj) = log(PxxFiltered+eps)*cos(jj.*([1:size(PxxFiltered,2)]' - 0.5).*(pi/nCepstralFeatures));
            end
            Cepstrum = Cepstrum(:,1:end-1);
            Cepstrum = reshape(Cepstrum,[size(Cepstrum,1)*size(Cepstrum,2),1])';

            
            obj.featuresVector = [obj.featuresVector,Cepstrum];
            
            
            
            function H = mel_spaced_filter_bank(TFRFreqs,NBanks)
                
                MelRange = freq_2_mel([TFRFreqs(1),TFRFreqs(end)]);
                Mels = linspace(MelRange(1),MelRange(2),NBanks+2);
                Freqs = mel_2_freq(Mels);
                location = zeros(size(Freqs));
                %Find nearest frequencies in TFR to target freqs
                for k = 1:length(Freqs)
                    [~, index] = min(abs(Freqs(k)-TFRFreqs));
                    location(k) = index;
                end
                locationUnique = unique(location);
                Freqs = TFRFreqs(locationUnique);
                NBanks = length(Freqs)-2;
                
                H = zeros(NBanks,length(TFRFreqs));
                for m = 2:NBanks+1
                    for k = 1:length(TFRFreqs)
                        H(m-1,k) = create_mel_filter_bank(k,locationUnique,m);
                    end
                end
                
                % figure
                % for m = 1:size(H,1)
                %     plot(TFRFreqs,H(m,:))
                %     hold on
                % end
                
                function freq = mel_2_freq(mel)
                    freq = 700.*(exp(mel./1125) - 1);
                end
                
                function mel = freq_2_mel(freq)
                    mel = 1125.*log(1 + freq./700);
                end
                
                function H=create_mel_filter_bank(k,f,m)
                    if k<f(m-1)
                        H = 0;
                    elseif (k>=f(m-1)&&k<=f(m))
                        H = (k-f(m-1))/(f(m)-f(m-1));
                    elseif (k>=f(m)&&k<=f(m+1))
                        H = (f(m+1)-k)/(f(m+1)-f(m));
                    elseif k>f(m+1)
                        H = 0;
                    end
                end
                
                
                
            end

            
            
        end
        
        %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Interbeat Features %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function obj = interbeat_features(obj)
            
            m_RR        = round(mean(diff(obj.assignedStates(:,1))));             % mean value of RR intervals
            sd_RR       = round(std(diff(obj.assignedStates(:,1))));              % standard deviation (SD) value of RR intervals
            mean_IntS1  = round(mean(obj.assignedStates(:,2)-obj.assignedStates(:,1)));            % mean value of S1 intervals
            sd_IntS1    = round(std(obj.assignedStates(:,2)-obj.assignedStates(:,1)));             % SD value of S1 intervals
            mean_IntS2  = round(mean(obj.assignedStates(:,4)-obj.assignedStates(:,3)));            % mean value of S2 intervals
            sd_IntS2    = round(std(obj.assignedStates(:,4)-obj.assignedStates(:,3)));             % SD value of S2 intervals
            mean_IntSys = round(mean(obj.assignedStates(:,3)-obj.assignedStates(:,2)));            % mean value of systole intervals
            sd_IntSys   = round(std(obj.assignedStates(:,3)-obj.assignedStates(:,2)));             % SD value of systole intervals
            mean_IntDia = round(mean(obj.assignedStates(2:end,1)-obj.assignedStates(1:end-1,4)));  % mean value of diastole intervals
            sd_IntDia   = round(std(obj.assignedStates(2:end,1)-obj.assignedStates(1:end-1,4)));   % SD value of diastole intervals
            
            for i=1:size(obj.assignedStates,1)-1
                R_SysRR(i)  = (obj.assignedStates(i,3)-obj.assignedStates(i,2))/(obj.assignedStates(i+1,1)-obj.assignedStates(i,1));
                R_DiaRR(i)  = (obj.assignedStates(i+1,1)-obj.assignedStates(i,4))/(obj.assignedStates(i+1,1)-obj.assignedStates(i,1));
                R_SysDia(i) = R_SysRR(i)/R_DiaRR(i);
                
                P_S1(i)     = sum(abs(obj.PCG(obj.assignedStates(i,1):obj.assignedStates(i,2))))/(obj.assignedStates(i,2)-obj.assignedStates(i,1));
                P_Sys(i)    = sum(abs(obj.PCG(obj.assignedStates(i,2):obj.assignedStates(i,3))))/(obj.assignedStates(i,3)-obj.assignedStates(i,2));
                P_S2(i)     = sum(abs(obj.PCG(obj.assignedStates(i,3):obj.assignedStates(i,4))))/(obj.assignedStates(i,4)-obj.assignedStates(i,3));
                P_Dia(i)    = sum(abs(obj.PCG(obj.assignedStates(i,4):obj.assignedStates(i+1,1))))/(obj.assignedStates(i+1,1)-obj.assignedStates(i,4));
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
            
            m_Ratio_SysRR   = mean(R_SysRR);
            sd_Ratio_SysRR  = std(R_SysRR);
            m_Ratio_DiaRR   = mean(R_DiaRR);
            sd_Ratio_DiaRR  = std(R_DiaRR);
            m_Ratio_SysDia  = mean(R_SysDia);
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
            
            interBeatFeatures = [m_RR, sd_RR, mean_IntS1, sd_IntS1, mean_IntS2, sd_IntS2,...
                mean_IntSys, sd_IntSys, mean_IntDia, sd_IntDia, m_Ratio_SysRR, sd_Ratio_SysRR,...
                m_Ratio_DiaRR, sd_Ratio_DiaRR, m_Ratio_SysDia, sd_Ratio_SysDia,m_Amp_SysS1,sd_Amp_SysS1,...
                m_Amp_DiaS2, sd_Amp_DiaS2];
            
            obj.featuresVector = [obj.featuresVector,interBeatFeatures];
            
        end
        
        
        
        %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Time domain complexity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        function obj = time_domain_complexity_features(obj,m,tau)
            
            middleBeat = round(size(obj.assignedStates,1)/2);
            PCGMiddleBeat = obj.PCG(obj.assignedStates(middleBeat,1):obj.assignedStates(middleBeat+1,1));
            PCGMiddleBeatPrime = diff(PCGMiddleBeat);
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
            
            
            zeroCrossingPoints = 0;
            for i = 1:(length(PCGMiddleBeat)-1)
                if PCGMiddleBeat(i)*PCGMiddleBeat(i+1) < 0
                    zeroCrossingPoints = zeroCrossingPoints+1;
                end
            end
            turningPoints = 0;
            for i = 1:(length(PCGMiddleBeatPrime)-1)
                if PCGMiddleBeatPrime(i)*PCGMiddleBeatPrime(i+1) < 0
                    turningPoints = turningPoints+1;
                end
            end
            
            simplicityVector = reshape(simplicityMatrix,[size( simplicityMatrix,1)*size( simplicityMatrix,2),1])';
            approximateEntropyVector = reshape(approximateEntropy,[size( approximateEntropy,1)*size( approximateEntropy,2),1])';
            
            obj.featuresVector = [obj.featuresVector,simplicityVector,approximateEntropyVector,zeroCrossingPoints,turningPoints];
            
        end
        
        %%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Spectral complexity %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = spectral_complexity_features(obj,nFrequencyBins,fractionOverlap,NS1,NSystole,NS2,NDiastole)
            
            
            [S1Segments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,1:2),NS1);
            [SystoleSegments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,2:3),NSystole);
            [S2Segments] = Features.get_timings_of_selected_segments(obj.assignedStates(:,3:4),NS2);
            [DiastoleSegments] = Features.get_timings_of_selected_segments([obj.assignedStates(1:end-1,4) obj.assignedStates(2:end,1)],NDiastole);
            
            [PxxSpecS,~] = Features.power_spectrum_MFCC(obj.PCG,S1Segments,SystoleSegments,S2Segments,DiastoleSegments,fractionOverlap,nFrequencyBins,obj.Fs);
            
            PxxSpecS = PxxSpecS./(max(max(PxxSpecS)));
            
            SpecS = -sum(PxxSpecS.*log(PxxSpecS),1);
            StandardDeviation = std(PxxSpecS,0,2)';
            Skewness = skewness(PxxSpecS,0,2)';
            Kurtosis = kurtosis(PxxSpecS,0,2)';
            
            obj.featuresVector = [obj.featuresVector, SpecS, StandardDeviation, Skewness, Kurtosis];
            
        end
        
    end
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    methods(Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function low_pass_filtered_signal = butterworth_low_pass_filter(original_signal,order,cutoff,sampling_frequency)
            %Get the butterworth filter coefficients
            [B_low,A_low] = butter(order,2*cutoff/sampling_frequency,'low');
            %Forward-backward filter the original signal using the butterworth
            %coefficients, ensuring zero phase distortion
            low_pass_filtered_signal = filtfilt(B_low,A_low,original_signal);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function high_pass_filtered_signal = butterworth_high_pass_filter(original_signal,order,cutoff,sampling_frequency)
            %Get the butterworth filter coefficients
            [B_high,A_high] = butter(order,2*cutoff/sampling_frequency,'high');
            %Forward-backward filter the original signal using the butterworth
            %coefficients, ensuring zero phase distortion
            high_pass_filtered_signal = filtfilt(B_high,A_high,original_signal);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [despiked_signal] = schmidt_spike_removal(original_signal, fs)
            
            %% Find the window size
            % (500 ms)
            windowsize = round(fs/2);
            
            %% Find any samples outside of a integer number of windows:
            trailingsamples = mod(length(original_signal), windowsize);
            
            %% Reshape the signal into a number of windows:
            sampleframes = reshape( original_signal(1:end-trailingsamples), windowsize, []);
            
            %% Find the MAAs:
            MAAs = max(abs(sampleframes));
            
            
            % While there are still samples greater than 3* the median value of the
            % MAAs, then remove those spikes:
            cont = 0;
            while( ~isempty(find((MAAs>median(MAAs)*3))) && cont < 1000)
                cont = cont + 1;
                %Find the window with the max MAA:
                [val window_num] = max(MAAs);
                if(numel(window_num)>1)
                    window_num = window_num(1);
                end
                
                %Find the postion of the spike within that window:
                [val spike_position] = max(abs(sampleframes(:,window_num)));
                
                if(numel(spike_position)>1)
                    spike_position = spike_position(1);
                end
                
                
                % Finding zero crossings (where there may not be actual 0 values, just a change from positive to negative):
                zero_crossings = [abs(diff(sign(sampleframes(:,window_num))))>1; 0];
                
                %Find the start of the spike, finding the last zero crossing before
                %spike position. If that is empty, take the start of the window:
                spike_start = max([1 find(zero_crossings(1:spike_position),1,'last')]);
                
                %Find the end of the spike, finding the first zero crossing after
                %spike position. If that is empty, take the end of the window:
                zero_crossings(1:spike_position) = 0;
                spike_end = min([(find(zero_crossings,1,'first')) windowsize]);
                
                %Set to Zero
                sampleframes(spike_start:spike_end,window_num) = 0.0001;
                
                %Recaclulate MAAs
                MAAs = max(abs(sampleframes));
            end
            if cont>999
                'This recording has had issues with spike removal'
            end
            
            despiked_signal = reshape(sampleframes, [],1);
            
            % Add the trailing samples back to the signal:
            despiked_signal = [despiked_signal; original_signal(length(despiked_signal)+1:end)];
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function normalizedTFR = per_frequency_normalization(TFR)
            normalizedTFR = zeros(size(TFR));
            for i = 1:size(TFR,2)
                tempScale = TFR(:,i);
                normalizedTFR(:,i) =  normalise_signal(tempScale);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [Segments] = get_timings_of_selected_segments(relevantStarts,N)
            Segments = cell(1,N);
            for i = 1:N
                Diff = (relevantStarts(:,2)) - relevantStarts(:,1);
                Segments{i} = int64([relevantStarts(:,1)+(Diff.*((i-1)./N)) relevantStarts(:,1)+(Diff.*((i)./N))-1]);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function FeatureOverBeat = get_feature_over_beat(Feature,S1Segments,SystoleSegments,S2Segments,DiastoleSegments)
            FeatureOverBeatsS1 = Features.mean_in_segment(Feature,S1Segments);
            FeatureOverBeatsSystole = Features.mean_in_segment(Feature,SystoleSegments);
            FeatureOverBeatsS2 = Features.mean_in_segment(Feature,S2Segments);
            FeatureOverBeatsDiastole = Features.mean_in_segment(Feature,DiastoleSegments);
            
            %NeedToRemoveFinalRow of S1/Systole/S2 as Diastole has 1 shorter
            FeatureOverBeatsS1 = FeatureOverBeatsS1(1:end-1,:);
            FeatureOverBeatsSystole = FeatureOverBeatsSystole(1:end-1,:);
            FeatureOverBeatsS2 = FeatureOverBeatsS2(1:end-1,:);
            
            FeatureOverBeats = horzcat(FeatureOverBeatsS1,FeatureOverBeatsSystole,FeatureOverBeatsS2,FeatureOverBeatsDiastole);
            
            % FeatureOverBeat = selectCorrelatedBeatsAndMean(FeatureOverBeats);
            FeatureOverBeat = FeatureOverBeats;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function FeatureOverBeat = mean_in_segment(Feature,Segments)
            FeatureOverBeat = zeros(size(Segments{1},1),length(Segments));
            for i = 1:length(Segments)
                for j = 1:size(Segments{1},1)
                    FeatureOverBeat(j,i) = mean(Feature(Segments{i}(j,1):Segments{i}(j,2),:),1);
                end
            end
            
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function meanOneBeat = select_correlated_beats_and_mean(ManyBeats,dimensionToAverageOver)
            if length(size(ManyBeats)) ~= 3 && length(size(ManyBeats)) ~= 2
                error('This method can only average over 2 or 3 dimensional feature vectors')
            end
            
            if length(size(ManyBeats)) == 3
                total = size(ManyBeats,1).*size(ManyBeats,2).*size(ManyBeats,3);
                %for a 3d vector concatenate the two dimensions you are not averaging
                %over
                eachBeat = zeros(size(ManyBeats,dimensionToAverageOver),total./size(ManyBeats,dimensionToAverageOver));
                for i = 1:size(ManyBeats,dimensionToAverageOver)
                    if dimensionToAverageOver == 1
                        eachBeat(i,:) =  reshape(ManyBeats(i,:,:),[1, total./size(ManyBeats,dimensionToAverageOver)]);
                    elseif dimensionToAverageOver == 2
                        eachBeat(i,:) =  reshape(ManyBeats(:,i,:),[1, total./size(ManyBeats,dimensionToAverageOver)]);
                    elseif dimensionToAverageOver == 3
                        eachBeat(i,:) =  reshape(ManyBeats(:,:,i),[1, total./size(ManyBeats,dimensionToAverageOver)]);
                    else
                        error('Dimension to average over must be 1,2,3')
                    end
                    
                end
                
            elseif length(size(ManyBeats)) == 2
                if dimensionToAverageOver == 2
                    eachBeat = ManyBeats';
                else
                    eachBeat = ManyBeats;
                end
            end
            
            
            Closest = pdist(eachBeat);
            Z = squareform(Closest);
            Z(Z == 0) = NaN;
            [minCol, iCol] = min(Z);
            [minValue, iRow] = min(minCol);
            FirstPoint = iRow;
            SecondPoint = iCol(iRow);
            extraPoint = [];
            for k = 1:length(Z)
                if k ~= FirstPoint || k~=SecondPoint
                    if (Z(FirstPoint,k) < 1.5*minValue && Z(SecondPoint,k) < 1.5*minValue)
                        extraPoint = [extraPoint k];
                    end
                end
            end
            HeartBeatsSelected = [extraPoint FirstPoint SecondPoint];
            
            eachBeat = eachBeat(HeartBeatsSelected,:);
            eachBeatMean = mean(eachBeat,1);
            
            if length(size(ManyBeats)) == 3
                if dimensionToAverageOver == 1
                    meanOneBeat = reshape(eachBeatMean,[size(ManyBeats,2),size(ManyBeats,3)]);
                elseif dimensionToAverageOver == 2
                    meanOneBeat = reshape(eachBeatMean,[size(ManyBeats,1),size(ManyBeats,3)]);
                elseif dimensionToAverageOver == 3
                    meanOneBeat = reshape(eachBeatMean,[size(ManyBeats,1),size(ManyBeats,2)]);
                else
                    error('Dimension to average over must be 1,2,3')
                end
                
                
            elseif length(size(ManyBeats)) == 2
                meanOneBeat = eachBeatMean;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Pxx,freqs] = power_spectrum_MFCC(signal,S1Segments,SystoleSegments,S2Segments,DiastoleSegments,fractionOverlap,nFrequencyBins,Fs)
            % A function designed to find the power of a signal in different time
            % bins given by S1Segments, SystoleSegments etc... tand different
            % frequency bins
            freqs = (Fs./2).*(linspace(0,1,nFrequencyBins+1));
            freqs = freqs(2:end);
            [S1Pxx] = pxx_from_single_segment(signal,S1Segments,fractionOverlap,freqs,Fs);
            [SystolePxx] = pxx_from_single_segment(signal,SystoleSegments,fractionOverlap,freqs,Fs);
            [S2Pxx] = pxx_from_single_segment(signal,S2Segments,fractionOverlap,freqs,Fs);
            [DiastolePxx] = pxx_from_single_segment(signal,DiastoleSegments,fractionOverlap,freqs,Fs);
            
            nTimeDomain = size(S1Pxx{1},2) + size(SystolePxx{1},2) + size(S2Pxx{1},2) + size(DiastolePxx{1},2);
            %     Pxx = zeros(nFrequencyBins,nTimeDomain);
            
            
            %             PxxAllBeats = zeros(size(S1Segments{1},1)-1,nFrequencyBins,nTimeDomain);
            for i = 1:nFrequencyBins
                FeatureOverBeats = horzcat(S1Pxx{i}(1:end-1,:),SystolePxx{i}(1:end-1,:),S2Pxx{i}(1:end-1,:),DiastolePxx{i});
                %                 PxxAllBeats(:,i,:) = FeatureOverBeats;
                FeatureOverBeat = Features.select_correlated_beats_and_mean(FeatureOverBeats,1);
                Pxx(i,:) = FeatureOverBeat;
            end
            %             Pxx = Features.select_correlated_beats_and_mean(PxxAllBeats,1);
            
            
            
            function pxxCell = pxx_from_single_segment(signal,segments,fractionOverlap,freqs,Fs)
                % A function to give the Pxx for a given segment (S1,Systole,S2 etc).
                % pxxCell has a length given by the number of frequency bins. Rumber of
                % columns is number of points in the segment (e.g. 7 time points in systole
                % is common and 3 in S1). Number of rows given by number of heart beats
                
                pxxCell = cell(length(freqs),1);
                for k = 1:length(freqs)
                    pxxCell{k} = zeros(size(segments{1},1),length(segments));
                end
                for j = 1:length(segments)
                    for ii = 1:size(segments{j},1)
                        startPointTemp = segments{j}(ii,1);
                        endPointTemp = segments{j}(ii,2);
                        startPoint = floor(startPointTemp - (endPointTemp-startPointTemp)*fractionOverlap);
                        endPoint = ceil(endPointTemp + (endPointTemp-startPointTemp)*fractionOverlap);
                        if startPoint < 1
                            startPoint = 1;
                        end
                        if endPoint > length(signal)
                            endPoint = length(signal);
                        end
                        x = signal(startPoint:endPoint);
                        %The NaN is in there as if you put a single number into
                        %periodogram (e.g 100) it calculate for 100 frequencies not at
                        %100 Hz which is what we want here
                        [pxxTemp] = periodogram(x,hamming(length(x)),freqs,Fs,'power');
                        for k = 1:length(freqs)
                            pxxCell{k}(ii,j) = abs(pxxTemp(k));
                        end
                    end
                end
            end
        end
    end
    
    
end
