function FeatureOverBeat = mean_in_segment(Feature,Segments)
FeatureOverBeat = zeros(size(Segments{1},1),length(Segments));
for i = 1:length(Segments)
    for j = 1:size(Segments{1},1)
        FeatureOverBeat(j,i) = mean(Feature(Segments{i}(j,1):Segments{i}(j,2),:),1);
    end
end

end