function [Pxx, freqs] = powerSpectrumMFCC(signal,S1Segments,SystoleSegments,S2Segments,DiastoleSegments,fractionOverlap,nFrequencyBins,Fs)
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
    Pxx = zeros(nFrequencyBins,nTimeDomain);
    for i = 1:nFrequencyBins
        FeatureOverBeats = horzcat(S1Pxx{i}(1:end-1,:),SystolePxx{i}(1:end-1,:),S2Pxx{i}(1:end-1,:),DiastolePxx{i});
        FeatureOverBeat = selectCorrelatedBeatsAndMean(FeatureOverBeats);
        Pxx(i,:) = FeatureOverBeat;
    end
end