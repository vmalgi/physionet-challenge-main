function FeatureOverBeat = get_feature_over_beat(Feature,S1Segments,SystoleSegments,S2Segments,DiastoleSegments)
FeatureOverBeatsS1 = Features.mean_in_segment(Feature,S1Segments);
FeatureOverBeatsSystole = Features.mean_in_segment(Feature,SystoleSegments);
FeatureOverBeatsS2 = Features.mean_in_segment(Feature,S2Segments);
FeatureOverBeatsDiastole = Features.mean_in_segment(Feature,DiastoleSegments);

%NeedToRemoveFinalRow of S1/Systole/S2 as Diastole has 1 shorter
FeatureOverBeatsS1 = FeatureOverBeatsS1(1:end-1,:);
FeatureOverBeatsSystole = FeatureOverBeatsSystole(1:end-1,:);
FeatureOverBeatsS2 = FeatureOverBeatsS2(1:end-1,:);

FeatureOverBeats = horzcat(FeatureOverBeatsS1,FeatureOverBeatsSystole,FeatureOverBeatsS2,FeatureOverBeatsDiastole);

% FeatureOverBeat = selectCorrelatedBeatsAndMean(FeatureOverBeats);
FeatureOverBeat = FeatureOverBeats;
end
