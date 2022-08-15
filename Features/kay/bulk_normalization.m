function normalizedTFR = bulk_normalization(TFR)
rowVectorTFR = reshape(TFR,1,[]);
normalizedRowVectorTFR = normalise_signal(rowVectorTFR);
normalizedTFR = reshape(normalizedRowVectorTFR,size(TFR));
end