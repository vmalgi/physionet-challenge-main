function [featuresVector] = getFeaturesVector(varargin)
%reorder features vector so that you have each feature over heartbeat
featuresVector = [];
for j = 1:nargin
    for i = 1:length(varargin{j})
        featuresVector = [featuresVector, varargin{j}{i}];
    end
end

end