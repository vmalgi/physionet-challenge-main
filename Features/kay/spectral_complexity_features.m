function featuresVector = spectral_complexity_features(PCG,Fs,A,nFrequencyBins,fractionOverlap,NS1,NSystole,NS2,NDiastole)
% spectral complexity
if nargin<4
    nFrequencyBins=5;
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

[PxxSpecS,~] = power_spectrum_MFCC(PCG,S1Segments,SystoleSegments,S2Segments,DiastoleSegments,fractionOverlap,nFrequencyBins,Fs);

PxxSpecS = PxxSpecS./(max(max(PxxSpecS)));

SpecS = -sum(PxxSpecS.*log(PxxSpecS),1);
StandardDeviation = std(PxxSpecS,0,2)';
Skewness = skewness(PxxSpecS,0,2)';
Kurtosis = kurtosis(PxxSpecS,0,2)';

featuresVector = [SpecS, StandardDeviation, Skewness, Kurtosis];

end
