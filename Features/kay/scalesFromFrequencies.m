function scales = scalesFromFrequencies(wavelet,freqrange,voicesPerOctave,fs)

% switch wavelet
%     case 'morl'
%         sigma = 5;
%         minscale = (sigma)./(2.*pi.*freqrange(2)*(1./fs));
%         maxscale = (sigma)./(2.*pi.*freqrange(1)*(1./fs));
%     case 'bump'
%         mu = 3;
%         sigma = 1.5;
%         minscale = (mu)./(2.*pi.*freqrange(2)*(1./fs));
%         maxscale = (mu)./(2.*pi.*freqrange(1)*(1./fs));
% end

cfreq = centfrq(wavelet);

minscale = cfreq./(freqrange(2).*(1./fs));
maxscale = cfreq./(freqrange(1).*(1./fs));

tempMinScale = floor(voicesPerOctave*log2(minscale));
tempMaxScale = ceil(voicesPerOctave*log2(maxscale));

scales = ((2^(1/voicesPerOctave)).^(tempMinScale:tempMaxScale));