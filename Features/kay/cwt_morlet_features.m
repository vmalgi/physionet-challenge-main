function featuresVector = cwt_morlet_features(PCG,Fs,A,voicesPerOctave,NS1,NSystole,NS2,NDiastole,freqRange)
% Continuous wavelet transform
if nargin<4
    voicesPerOctave=2.4;
    NS1=3;
    NSystole=7;
    NS2=3;
    NDiastole=7;
    freqRange=[30,450];
end

wavelet_name = 'morl';
scales = scales_from_frequencies(wavelet_name,freqRange,voicesPerOctave,Fs);
if(length(PCG)< Fs*1.025)
    PCG = [PCG; zeros(round(0.025*Fs),1)];
end
[cd] = cwt(PCG,scales,wavelet_name);
cd = cd';
TFR = abs(cd);
normalizedTFR = per_frequency_normalization(TFR);

NTot = NS1 + NSystole + NS2 + NDiastole;

[S1Segments] = get_timings_of_selected_segments(A(:,1:2),NS1);
[SystoleSegments] = get_timings_of_selected_segments(A(:,2:3),NSystole);
[S2Segments] = get_timings_of_selected_segments(A(:,3:4),NS2);
[DiastoleSegments] = get_timings_of_selected_segments([A(1:end-1,4) A(2:end,1)],NDiastole);

meanCWTAllBeats = zeros(size(A,1)-1,NTot,size(normalizedTFR,2));
meanCWTMatrix = zeros(NTot,size(normalizedTFR,2));
for i = 1:size(normalizedTFR,2)
    meanCWTAllBeats(:,:,i) = get_feature_over_beat(normalizedTFR(:,i),S1Segments,...
        SystoleSegments,S2Segments,DiastoleSegments);
    
    meanCWTMatrix(:,i) = select_correlated_beats_and_mean(meanCWTAllBeats(:,:,i),1);
end

meanCWT = reshape(meanCWTMatrix,[1,size(meanCWTMatrix,1)*size(meanCWTMatrix,2)]);

featuresVector = meanCWT;

    function scales = scales_from_frequencies(wavelet,freqrange,voicesPerOctave,fs)
        
        cfreq = centfrq(wavelet);
        
        minscale = cfreq./(freqrange(2).*(1./fs));
        maxscale = cfreq./(freqrange(1).*(1./fs));
        
        tempMinScale = floor(voicesPerOctave*log2(minscale));
        tempMaxScale = ceil(voicesPerOctave*log2(maxscale));
        
        scales = ((2^(1/voicesPerOctave)).^(tempMinScale:tempMaxScale));
        
    end
end