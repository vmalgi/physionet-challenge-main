function [Segments] = get_timings_of_selected_segments(relevantStarts,N)

Segments = cell(1,N);
for i = 1:N
    Diff = (relevantStarts(:,2)) - relevantStarts(:,1);
    Segments{i} = int64([relevantStarts(:,1)+(Diff.*((i-1)./N)) relevantStarts(:,1)+(Diff.*((i)./N))-1]); 
end

