function H = melSpacedFilterBank(TFR,TFRFreqs,NBanks)

MelRange = freq2mel([TFRFreqs(1),TFRFreqs(end)]); 
Mels = linspace(MelRange(1),MelRange(2),NBanks+2);
Freqs = mel2freq(Mels);
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
        H(m-1,k) = createMelFilterBank(k,locationUnique,m);
    end
end

% figure
% for m = 1:size(H,1)
%     plot(TFRFreqs,H(m,:))
%     hold on
% end



end

