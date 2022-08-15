function pxxCell = pxx_from_single_segment(signal,segments,fractionOverlap,freqs,Fs)
% A function to give the Pxx for a given segment (S1,Systole,S2 etc).
% pxxCell has a length given by the number of frequency bins. Rumber of
% columns is number of points in the segment (e.g. 7 time points in systole
% is common and 3 in S1). Number of rows given by number of heart beats

pxxCell = cell(length(freqs),1);
for k = 1:length(freqs)
    pxxCell{k} = zeros(length(segments{1}),length(segments));
end
for j = 1:length(segments)
    for i = 1:length(segments{j})
        startPointTemp = segments{j}(i,1);
        endPointTemp = segments{j}(i,2);
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
            pxxCell{k}(i,j) = abs(pxxTemp(k));
        end
    end
end


end