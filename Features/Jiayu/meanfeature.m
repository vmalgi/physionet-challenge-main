%% obtain the mean vector  
% special for our feature vector 1+2+1+4,
% respectively for S1, sys, S2, dia
function feature_mean = meanfeature(feature,timeN,freqN)
L = size(feature,2);
cycle = 8*timeN;
num = floor(L/cycle);
count_table = [1 2 1 4]*timeN;
feature_mean = [];
for i = 1:4
    featuretemp = [];
    for k = 1:num
        start = i+(k-1)*cycle;
        temp = reshape(feature(:,start:start+count_table(i)-1),1,freqN*count_table(i));
        featuretemp = [featuretemp; temp];
    end
    featuretemp = mean(featuretemp);
    feature_mean = [feature_mean featuretemp];
end


