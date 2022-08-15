function F = extractFeaturesFromPCG(ss,PCG)
%
% This function calculate features based on the assigned_states by running "runSpringerSegmentationAlgorithm.m" function
%
% INPUTS:
% assigned_states: the array of state values assigned to the sound recording.
% PCG: resampled sound recording at 1000 Hz.
%
% OUTPUTS:
% features: matrix with size #cycles x #features (36 features per segment)
%
%% Start with first complete cycle
n=length(ss);
fnzs=min(find(abs(PCG)>1E-5)); % skip leading "zeros" signal
ss=ss(fnzs+1:end);
PCG=PCG(fnzs+1:end);
indx = find(abs(diff(ss))>0)+1;
indx=[1;indx;n-fnzs];
% start at first S1 cycle
f1=min(find(ss==1)); % first complete cycle
fcy1=min(find(indx==f1));
indx=indx(fcy1:end);
Ncy=floor(length(indx)/4);
F=[];
%% Feature calculation only for complete cycles
for i=0:Ncy-2
    S1=indx(4*i+2)-indx(4*i+1);
    S2=indx(4*i+3)-indx(4*i+2);
    S3=indx(4*i+4)-indx(4*i+3);
    S4=indx(4*i+5)-indx(4*i+4);
    %
    scale=1;%max(mad(PCG(indx(4*i+1):indx(4*i+5))),1e-10);
    PCG1n=(PCG(indx(4*i+1):indx(4*i+2)))/scale;
    PCG2n=(PCG(indx(4*i+2):indx(4*i+3)))/scale;
    PCG3n=(PCG(indx(4*i+3):indx(4*i+4)))/scale;
    PCG4n=(PCG(indx(4*i+4):indx(4*i+5)))/scale;
    %
   if S1>=256
        MLCCS1=melcepst(PCG1n,1000,'EdDNny');
    else
        MLCCS1=melcepst(wextend(1,'per',PCG1n,256-S1,'r'),1000,'EdDNny');
    end
    if S2>=256
        MLCCS2=melcepst(PCG2n,1000,'EdDNny');
    else
        MLCCS2=melcepst(wextend(1,'per',PCG2n,256-S2,'r'),1000,'EdDNny');
    end
    if S3>=256
        MLCCS3=melcepst(PCG3n,1000,'EdDNny');
    else
        MLCCS3=melcepst(wextend(1,'per',PCG3n,256-S3,'r'),1000,'EdDNny');
    end
    if S4>=256
        MLCCS4=melcepst(PCG4n,1000,'EdDNny');
    else
        MLCCS4=melcepst(wextend(1,'per',PCG4n,256-S4,'r'),1000,'EdDNny');
    end
    %
    % Wavelets mad
    %
    N=6;
    NameWT='sym8';
if length(PCG1n)<64
        PCG1n=wextend(1,'per',PCG1n,64-S1,'r');
    end
    if length(PCG2n)<64
        PCG2n=wextend(1,'per',PCG2n,64-S2,'r');
    end
    if length(PCG3n)<64
        PCG3n=wextend(1,'per',PCG3n,64-S3,'r');
    end
    if length(PCG4n)<64
        PCG4n=wextend(1,'per',PCG4n,64-S4,'r');
    end
    WLTS1=modwpt(PCG1n,NameWT,N);
    WLTS2=modwpt(PCG2n,NameWT,N);
    WLTS3=modwpt(PCG3n,NameWT,N);
    WLTS4=modwpt(PCG4n,NameWT,N);
    WLTS1d=mad(WLTS1');
    WLTS2d=mad(WLTS2');
    WLTS3d=mad(WLTS3');
    WLTS4d=mad(WLTS4');
%
%
    MLCC1d=mad(MLCCS1);
    MLCC2d=mad(MLCCS2);
    MLCC3d=mad(MLCCS3);
    MLCC4d=mad(MLCCS4);
    %
    MLCC1m=mean(MLCCS1);
    MLCC2m=mean(MLCCS2);
    MLCC3m=mean(MLCCS3);
    MLCC4m=mean(MLCCS4);
%
    A=[MLCC1d MLCC1m WLTS1d MLCC2d MLCC2m WLTS2d MLCC3d MLCC3m WLTS3d MLCC4d MLCC4m WLTS4d];
    if ~any(isnan(A)) %do not include cycles with segments of "zero" signal
        F(:,i+1)=A;
    end
end
