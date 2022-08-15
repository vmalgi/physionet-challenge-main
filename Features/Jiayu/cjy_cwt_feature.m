%% CWT Feature Extraction
% we compute the following feature combinations
% isfilter   1 filter with 25-400Hz bandpass filter
%            0 without filter
% ismean     1 feature = sum(energy)/length(time)
%            0 feature = sum(energy)
% timeN      the number of time intervals S1:sys:S2:dia = 1:2:1:4
% feqN       the number of frequency intervals
function feature = cjy_cwt_feature(state,signal)

% paramters
Level = 32;
timeN = 2;
freqN = 4;
ismean = 0;
isfilter = 0;
energyind = [];
% go loop
if(isfilter)
    signal = butterworth_low_pass_filter(signal,2,400,Fs, false);
    signal = butterworth_high_pass_filter(signal,2,25,Fs);
end
ccfs=cwt(signal,1:Level,'db6');
idseg = find(diff(state) ~= 0);

% starts from the second complete cycle
p = find(state(idseg) == 1,2); % first 2
start = p(2)-1;
for j = start:length(idseg)-1 % reject the last segment
    energytemp = cwtSubfeature(ccfs(:,idseg(j)+1:idseg(j+1)),timeN,freqN,Level,state(idseg(j)+1),ismean);
    energyind = [energyind energytemp];
    %fprintf('the size of energyind is %d,the seg index is %d\n',size(energyind,2),statetemp(idseg(j)+1));
end
energy = meanfeature(energyind,timeN,freqN);
cwt_features = energy;

t_20_features = extractFeaturesFromHsIntervals(state,signal);
feature = [cwt_features t_20_features];