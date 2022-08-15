function featuresVector = mfcc_features(PCG,Fs,A,nCepstralFeatures,nFilterBanks,nFrequencyBins,fractionOverlap,NS1,NSystole,NS2,NDiastole)
if nargin<4
    nCepstralFeatures=21;
    nFilterBanks=20;
    nFrequencyBins=40;
    fractionOverlap=0.1;
    NS1=3;
    NSystole=7;
    NS2=3;
    NDiastole=7;
end


[S1Segments] = get_timings_of_selected_segments(A(:,1:2),NS1);
[SystoleSegments] = get_timings_of_selected_segments(A(:,2:3),NSystole);
[S2Segments] = get_timings_of_selected_segments(A(:,3:4),NS2);
[DiastoleSegments] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastole);

[Pxx,f] = power_spectrum_MFCC(PCG,S1Segments,SystoleSegments,S2Segments,DiastoleSegments,fractionOverlap,nFrequencyBins,Fs);

H = mel_spaced_filter_bank(f,nFilterBanks);


PxxFiltered = (H*Pxx)';

Cepstrum = zeros(size(PxxFiltered,1),nCepstralFeatures);
for jj = 1:nCepstralFeatures
    Cepstrum(:,jj) = log(PxxFiltered+eps)*cos(jj.*([1:size(PxxFiltered,2)]' - 0.5).*(pi/nCepstralFeatures));
end
Cepstrum = Cepstrum(:,1:end-1);
Cepstrum = reshape(Cepstrum,[size(Cepstrum,1)*size(Cepstrum,2),1])';


featuresVector = Cepstrum; 

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
