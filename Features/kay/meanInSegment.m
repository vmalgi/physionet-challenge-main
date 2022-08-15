
function FeatureOverBeat = meanInSegment(Feature,Segments)
FeatureOverBeat = zeros(length(Segments{1}),length(Segments));
for i = 1:length(Segments)
    for j = 1:length(Segments{i})
        FeatureOverBeat(j,i) = mean(Feature(Segments{i}(j,1):Segments{i}(j,2),:),1);
    end
end

end