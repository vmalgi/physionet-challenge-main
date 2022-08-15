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
