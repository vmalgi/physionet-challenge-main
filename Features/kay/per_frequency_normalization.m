function normalizedTFR = per_frequency_normalization(TFR)

normalizedTFR = zeros(size(TFR));
for i = 1:size(TFR,2)
    tempScale = TFR(:,i);
    normalizedTFR(:,i) =  normalise_signal(tempScale);
end

end
