%% My_mapminmax
% The built-in mapminmax maps data into [-1,1]
% which is incompatible with sigmoid.
% Therefore, this mapminmax maps data into [0,1]
function [dataset,maxref,range] = my_mapminmax(dataset)

maxref = max(dataset);
minref = min(dataset);

range = maxref-minref;
dataset = bsxfun(@minus,dataset,minref);
dataset = bsxfun(@rdivide,dataset,range);