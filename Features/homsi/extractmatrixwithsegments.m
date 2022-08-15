function A = extractmatrixwithsegments(assigned_states)
%
% This function calculate 20 features based on the assigned_states by running "runSpringerSegmentationAlgorithm.m" function
%
% INPUTS:
% assigned_states: the array of state values assigned to the sound recording.
% PCG: resampled sound recording with 1000 Hz.
%
% OUTPUTS:
% features: the obtained 20 features for the current sound recording
%
%
% Written by: Chengyu Liu, January 22 2016
%             chengyu.liu@emory.efdu
%
% Last modified by:
%
%
% $$$$$$ IMPORTANT
% Please note: the calculated 20 features are only some pilot features, some features maybe
% helpful for classifying normal/abnormal heart sounds, some maybe
% not. You need re-construct the features for a more accurate classification.


%% We just assume that the assigned_states cover at least 2 whole heart beat cycle
indx = find(abs(diff(assigned_states))>0); % find the locations with changed states

if assigned_states(1)>0   % for some recordings, there are state zeros at the beginning of assigned_states
    switch assigned_states(1)
        case 4
            K=1;
        case 3
            K=2;
        case 2
            K=3;
        case 1
            K=4;
    end
else
    switch assigned_states(indx(1)+1)
        case 4
            K=1;
        case 3
            K=2;
        case 2
            K=3;
        case 1
            K=0;
    end
    K=K+1;
end

indx2                = indx(K:end);
rem                  = mod(length(indx2),4);
indx2(end-rem+1:end) = [];
A                    = reshape(indx2,4,length(indx2)/4)'; % A is N*4 matrix, the 4 columns save the beginnings of S1, systole, S2 and diastole in the same heart cycle respectively
