%div_segf = amie_seg(Xabs, N, S, fs);

function [div_segf] = amie_seg_test(E, N, S, fs)
%
%
% Objetive:  this function implements the respiratory stage segmenter defined 
%            in reference "Automatic Multi-level In-Exhale Segmentation and Enhanced 
%            Generalized S-Transform for whezing detection". Specifically, this function
%            allows to segment the mixture spectrogram X into inspiratory and expiratory stages.
%
%
% Input: 
% - Xabs:    mixture spectrogram
% - N:       hamming window sample length
% - S:       overlap (between 0 and 1)
% - fs:      sample rate (Hz)
%
% Output: 
% - div_segf: vector that indicates the separation frames between segments.
%
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------

%% Initial thresholds (Th1, Th2)
Th1 = max(E)/4;
Th2 = mean(E)/2;
%% Adaptive threshold
Th_ad = 0.5;
Th_add = zeros(1,100); 
Th_add(1) = Th_ad;
cont = 2;
while Th_ad<Th2
    Th_ad = Th_ad + abs(Th1-Th2);
    Th_add(cont) = Th_ad;
    cont = cont + 1;
    if cont>100
        break
    end
end
Numb = length(find(Th_add>0));
Th_ad = Th_add(1:Numb);
%% Annotation (Silent Frame, Potential-breath Frame and Breath Frame)
Flag_S = zeros(Numb,length(E));
Flag_P = zeros(Numb,length(E));
Flag_B = zeros(Numb,length(E));
% Flagging frames
for i = 1:Numb
    %Silent frames
    Flag_S(i,:) = E<Th_ad(i);
    % Potential breath frames
    Flag_P(i,:) = E>Th_ad(i) & E<Th1;
    % Breath frames
    Flag_B(i,:) = E>Th1;
end
% Short cut modified method
% Finds potential breath and breath frames 
if length(Th_ad)>=3
decision = E>Th_ad(3);
decision =decision';
else
decision = E>Th_ad(2);
decision =decision';
end
%% Identification of the initial and final frame of each segment
POI = diff(decision);
start_points = (find(POI==1));
end_points = (find(POI==-1)); 
% Ensuring start is at the beginning and end is at the end
if (end_points(1) < start_points(1))
    start_points = [1; start_points];
end
if (end_points(end) < start_points(end))
    end_points = [end_points; length(decision)+1];
end
locations = [start_points end_points];
% minimum length of segments is 31.25ms
thres=31.25/1000*fs; 
locations = locations((locations(:,2)-locations(:,1))>thres,:);
%% Segmentation
div_seg=ones(1,size(locations,1)-1); 
for k = 1:size(locations,1)-1
   div_seg(k) =  locations(k,2)+round((locations(k+1,1)-locations(k,2))/2);
end
%%  Annotation 
div_segf=[1,div_seg,size(E,2)]; 
div_s = div_segf/fs;
end


