function features = get_heart_segmentation_features(assigned_states, PCG)
%% Paper Information
% Written by: Chengyu Liu, January 22 2016 chengyu.liu@emory.edu
%% Purpose
% Extract features based on assigned states
%% Input
% assigned_states: the array of state values assigned to the sound recording.
% PCG: resampled sound recording with 1000 Hz.
%% Outputs
% features
features=table; 

% find the locations with changed states
indx = find(abs(diff(assigned_states))>0);

% for some recordings, there are state zeros at the beginning of assigned_states
if assigned_states(1)>0
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

indx2 = indx(K:end);
rem = mod(length(indx2),4);
indx2(end-rem+1:end) = [];
% A is N*4 matrix, the 4 columns save the beginnings of S1, systole, S2 and diastole in the same heart cycle respectively
A = reshape(indx2,4,length(indx2)/4)';

%% Feature calculation
%Abdollahpur
[signal_s1,signal_s2, signal_systole,signal_diastole]=seperate_states_2(assigned_states,PCG,false);
a=mofify_A(A,signal_systole,signal_diastole,signal_s1,signal_s2);
if a(end)==0
    a(end)=length(PCG); 
end
% extracting  11 wtentropy features
%features.entropy_matrix=mean_beat_wt_entropy(PCG,a);
% extracting 9 wt features 
features.sum_matrix_mean=absolute_sum(PCG,a);
% extracting 9 shannon featueres 
features.shannon_matrix=shannon_en(PCG,a);
% calculating mfcc features
%[~,mean_mfcc_systole ,mean_mfcc_diastole]=calculating_mfcc(signal_s1 ,signal_s2 ,signal_systole, signal_diastole,a);
%features.feature_taki_mfcc=max( mean_mfcc_systole(4) ,mean_mfcc_diastole(4) );
%Jiayu
features.cjy_cwt_feature = cjy_cwt_feature_modified(assigned_states,PCG);

% Kay
Fs=1000;
%features.spectral_complexity_features = spectral_complexity_features(PCG,Fs,A);
features.mfcc_features = mfcc_features(PCG,Fs,A);
features.cwt_morlet_features = cwt_morlet_features(PCG,Fs,A);

% Bobillo
F = extractFeaturesFromPCG_modified(assigned_states,PCG);
% Pad with zeros to complete 172 cycles, if shorter.
st=size(F);
if st(2)<172
    Fp=F;
    Fp(:,st(2)+1:172,:)=0;
else
    % cut tail if longer
    Fp=F(:,1:172,:);
end
% duplicate to avoid recording for singleton 4-way tensor
Ft=cat(3,Fp,Fp);
% Reshape into 4-way tensor
for i=1:4
    T4(1:142,i,:,:)=Ft(142*(i-1)+1:142*i,:,:);
end
% Fill trailing zeros with copies of leading non-zeros
T4=fillTensor4(T4);
load('Params.mat', 'indx','U')
% Run classifier
f=grabFeatures4(T4,U,indx);
features.bobillo_features=f(1,:);


%% Time intervals
% % mean value of RR intervals (in samples)
% features.m_RR        = round(mean(diff(A(:,1))));
% % standard deviation (SD) value of RR intervals
% features.sd_RR       = round(std(diff(A(:,1))));
% % mean value of S1 intervals
% features.mean_IntS1  = round(mean(A(:,2)-A(:,1)));
% % SD value of S1 intervals
% features.sd_IntS1    = round(std(A(:,2)-A(:,1)));
% % mean value of S2 intervals
% features.mean_IntS2  = round(mean(A(:,4)-A(:,3)));
% % SD value of S2 intervals
% features.sd_IntS2    = round(std(A(:,4)-A(:,3)));
% % mean value of systole intervals
% features.mean_IntSys = round(mean(A(:,3)-A(:,2)));
% % SD value of systole intervals
% features.sd_IntSys   = round(std(A(:,3)-A(:,2)));
% % mean value of diastole intervals
% features.mean_IntDia = round(mean(A(2:end,1)-A(1:end-1,4)));
% % SD value of diastole intervals
% features.sd_IntDia   = round(std(A(2:end,1)-A(1:end-1,4)));

% homsi
N_LEVELS = 5;
WAVELET_NAME = 'db4';
[C,L] = wavedec(PCG, N_LEVELS, WAVELET_NAME);

rD=cell(1,N_LEVELS);
rA=cell(1,N_LEVELS);
for iLevel = N_LEVELS:-1:1
    rD{iLevel} = wrcoef('d', C, L, WAVELET_NAME, iLevel);
    rA{iLevel} = wrcoef('a', C, L, WAVELET_NAME, iLevel);
end

% initalisation
col=size(A,1)-1;
R_SysRR= zeros(1,col);
R_DiaRR= zeros(1,col);
R_SysDia= zeros(1,col);
P_S1= zeros(1,col);
P_Sys= zeros(1,col);
P_S2= zeros(1,col);
P_Dia= zeros(1,col);
P_SysS1= zeros(1,col);
P_DiaS2= zeros(1,col);
SK_S1= zeros(1,col);
% SK_Sys= zeros(1,col);
% SK_S2= zeros(1,col);
% SK_Dia= zeros(1,col);
KU_S1= zeros(1,col);
% KU_Sys= zeros(1,col);
KU_S2= zeros(1,col);
% KU_Dia= zeros(1,col);

%zcr_S1= zeros(1,col);
%zcr_SYS= zeros(1,col);
%zcr_S2= zeros(1,col);
%zcr_DIA= zeros(1,col);
%T_S1= zeros(1,col);
%T_Sys= zeros(1,col);
%T_S2= zeros(1,col);
%T_Dia= zeros(1,col);
%Max_S1= zeros(1,col);
Max_Sys= zeros(1,col);
%Max_S2= zeros(1,col);
%Max_DIA= zeros(1,col);
%mean_S1= zeros(1,col);
%mean_SYS= zeros(1,col);
%mean_S2= zeros(1,col);
%mean_DIA= zeros(1,col);
%rms_S1= zeros(1,col);
rms_SYS= zeros(1,col);
%rms_S2= zeros(1,col);
rms_DIA= zeros(1,col);
%ptd_S1= zeros(1,col);
ptd_SYS= zeros(1,col);
%ptd_S2= zeros(1,col);
ptd_DIA= zeros(1,col);
%pfd_S1= zeros(1,col);
pfd_SYS= zeros(1,col);
%pfd_S2= zeros(1,col);
pfd_DIA= zeros(1,col);
% saen_S1= zeros(1,col);
% saen_SYS= zeros(1,col);
% saen_S2= zeros(1,col);
% saen_DIA= zeros(1,col);
pks_S1max= zeros(1,col);
index1= zeros(1,col);
pks_SYSmax= zeros(1,col);
index2= zeros(1,col);
pks_S2max= zeros(1,col);
index3= zeros(1,col);
pks_DIAmax= zeros(1,col);
index4= zeros(1,col);
% BWs_S1= zeros(1,col);
BWs_SYS= zeros(1,col);
% BWs_S2= zeros(1,col);
% BWs_DIA= zeros(1,col);
% Qs1= zeros(1,col);
Qsys= zeros(1,col);
% Qs2= zeros(1,col);
% Qdia= zeros(1,col);
Shannon_S1= zeros(1,col);
Shannon_SYS= zeros(1,col);
% Shannon_S2= zeros(1,col);
Shannon_DIA= zeros(1,col);
% var_S1= zeros(1,col);
var_SYS= zeros(1,col);
% var_S2= zeros(1,col);
var_DIA= zeros(1,col);
% Shannon_A5_s1= zeros(1,col);
% Shannon_A5_sys= zeros(1,col);
% Shannon_A5_s2= zeros(1,col);
% Shannon_A5_dia= zeros(1,col);
% Shannon_D5_s1= zeros(1,col);
% Shannon_D5_sys= zeros(1,col);
% Shannon_D5_s2= zeros(1,col);
% Shannon_D5_dia= zeros(1,col);
% Shannon_D4_s1= zeros(1,col);
Shannon_D4_sys= zeros(1,col);
Shannon_D4_s2= zeros(1,col);
% Shannon_D4_dia= zeros(1,col);
% Shannon_D3_s1= zeros(1,col);
Shannon_D3_sys= zeros(1,col);
% Shannon_D3_s2= zeros(1,col);
% Shannon_D3_dia= zeros(1,col);
Shannon_D2_s1= zeros(1,col);
Shannon_D2_sys= zeros(1,col);
Shannon_D2_s2= zeros(1,col);
Shannon_D2_dia= zeros(1,col);
% Shannon_D1_s1= zeros(1,col);
% Shannon_D1_sys= zeros(1,col);
% Shannon_D1_s2= zeros(1,col);
% Shannon_D1_dia= zeros(1,col);

%detailed_segmentation_features_S1=table; 
detailed_segmentation_features_Sys=table; 
detailed_segmentation_features_S2=table; 
detailed_segmentation_features_Dia=table; 
for i=1:size(A,1)-1
    R_SysRR(i)  = (A(i,3)-A(i,2))/(A(i+1,1)-A(i,1))*100;
    R_DiaRR(i)  = (A(i+1,1)-A(i,4))/(A(i+1,1)-A(i,1))*100;
    R_SysDia(i) = R_SysRR(i)/R_DiaRR(i)*100;
    
    %% fix later
    %detailed_segmentation_features_S1(i,:)  = get_detailed_segmentation_features(PCG(A(i,1):A(i,2)),Fs);
    detailed_segmentation_features_Sys(i,:) = get_detailed_segmentation_features_sys(PCG(A(i,2):A(i,3)),Fs);
    detailed_segmentation_features_S2(i,:)  = get_detailed_segmentation_features_s2(PCG(A(i,3):A(i,4)),Fs);
    detailed_segmentation_features_Dia(i,:) = get_detailed_segmentation_features_dia(PCG(A(i,4):A(i+1,1)),Fs);
    
    %skewness (potes)
    SK_S1(i)  = skewness(PCG(A(i,1):A(i,2)));
%     SK_Sys(i) = skewness(PCG(A(i,2):A(i,3)));
%     SK_S2(i)  = skewness(PCG(A(i,3):A(i,4)));
%     SK_Dia(i) = skewness(PCG(A(i,4):A(i+1,1)));
    
    % kurtosis (potes)
    KU_S1(i)  = kurtosis(PCG(A(i,1):A(i,2)));
%     KU_Sys(i) = kurtosis(PCG(A(i,2):A(i,3)));
    KU_S2(i)  = kurtosis(PCG(A(i,3):A(i,4)));
%     KU_Dia(i) = kurtosis(PCG(A(i,4):A(i+1,1)));
    
    P_S1(i)     = sum(abs(PCG(A(i,1):A(i,2))))/(A(i,2)-A(i,1));
    P_Sys(i)    = sum(abs(PCG(A(i,2):A(i,3))))/(A(i,3)-A(i,2));
    P_S2(i)     = sum(abs(PCG(A(i,3):A(i,4))))/(A(i,4)-A(i,3));
    P_Dia(i)    = sum(abs(PCG(A(i,4):A(i+1,1))))/(A(i+1,1)-A(i,4));
    if P_S1(i)>0
        P_SysS1(i) = P_Sys(i)/P_S1(i)*100;
    else
        P_SysS1(i) = 0;
    end
    if P_S2(i)>0
        P_DiaS2(i) = P_Dia(i)/P_S2(i)*100;
    else
        P_DiaS2(i) = 0;
    end
    
    %% Homsi
    %Signals of S1, Sys, S2 & Dia
    S1=PCG(A(i,1):A(i,2));      L_S1=size(S1);      LS1=L_S1(1,1);
    SYS=PCG(A(i,2):A(i,3));     L_SYS=size(SYS);    LSYS=L_SYS(1,1);
    S2=PCG(A(i,3):A(i,4));      L_S2=size(S2);      LS2=L_S2(1,1);
    DIA=PCG(A(i,4):A(i+1,1));   L_DIA=size(DIA);    LDIA=L_DIA(1,1);
    % Zero Crossing
%     zcr_S1(i) = sum(diff(S1)>0)/length(S1);
%     zcr_SYS(i) = sum(diff(SYS)>0)/length(SYS);
%     zcr_S2(i) = sum(diff(S2)>0)/length(S2);
%     zcr_DIA(i) = sum(diff(DIA)>0)/length(DIA);
    % Duration of S1, Sys, S2 & Dia
%     T_S1(i)=A(i,2)-A(i,1);
%     T_Sys(i)=A(i,3)-A(i,2);
%     T_S2(i)=A(i,4)-A(i,3);
%     T_Dia(i)=A(i+1,1)-A(i,4);
    % Maximum
%     Max_S1(i)=max(S1);
    Max_Sys(i)=max(SYS);
%     Max_S2(i)=max(S2);
%     Max_DIA(i)=max(DIA);
    % Mean
%     mean_S1(i)=mean(S1);
%     mean_SYS(i)=mean(SYS);
%     mean_S2(i)=mean(S2);
%     mean_DIA(i)=mean(DIA);
    % RMS
%     rms_S1(i)=rms(S1);
    rms_SYS(i)=rms(SYS);
%     rms_S2(i)=rms(S2);
    rms_DIA(i)=rms(DIA);
    % Total Power in Time Domain
%     ptd_S1(i)=(norm(S1)^2)/LS1;
    ptd_SYS(i)=(norm(SYS)^2)/LSYS;
%     ptd_S2(i)=(norm(S2)^2)/LS2;
    ptd_DIA(i)=(norm(DIA)^2)/LDIA;
    % Total Power in Frequency Domain
%     fft_S1=fft(S1);
    fft_SYS=fft(SYS);
%     fft_S2=fft(S2);
    fft_DIA=fft(DIA);
%     pfd_S1(i)=sum(fft_S1.*conj(fft_S1))/(LS1^2);
    pfd_SYS(i)=sum(fft_SYS.*conj(fft_SYS))/(LSYS^2);
%     pfd_S2(i)=sum(fft_S2.*conj(fft_S2))/(LS2^2);
    pfd_DIA(i)=sum(fft_DIA.*conj(fft_DIA))/(LDIA^2);
    % Sample Entropy
%     saen_S1(i)=SampEn(2, 0.2*std(S1), S1, 1);
%     saen_SYS(i)=SampEn(2, 0.2*std(SYS), SYS, 1);
%     saen_S2(i)=SampEn(2, 0.2*std(S2), S2, 1);
%     saen_DIA(i)=SampEn(2, 0.2*std(DIA), DIA, 1);
    
    % FRECUENCY DOMAIN
    % Next power of 2 from length of S1
    NFFT1 = 2^nextpow2(LS1);
    fft_S1_NFFT = fft(S1,NFFT1)/LS1;
    % Next power of 2 from length of SYS
    NFFT2 = 2^nextpow2(LSYS);
    fft_SYS_NFFT = fft(SYS,NFFT2)/LSYS;
    % Next power of 2 from length of SYS
    NFFT3 = 2^nextpow2(LS2);
    fft_S2_NFFT = fft(S2,NFFT3)/LS2;
    
    NFFT4 = 2^nextpow2(LDIA); % Next power of 2 from length of SYS
    fft_DIA_NFFT = fft(DIA,NFFT4)/LDIA;
    
    % conditional
    if i==1
        [pks_S1,~] = findpeaks(2*abs(fft_S1_NFFT(1:NFFT1/2+1)));
        [pks_SYS,~] = findpeaks(2*abs(fft_SYS_NFFT(1:NFFT2/2+1)));
        [pks_S2,~] = findpeaks(2*abs(fft_S2_NFFT(1:NFFT3/2+1)));
        [pks_DIA,~] = findpeaks(2*abs(fft_DIA_NFFT(1:NFFT4/2+1)));
        
        if isempty(pks_S1)
            pks_S1max(i)=0;
            
        else
            [pks_S1max(i), index1(i)] = max(pks_S1);
        end
        
        if isempty(pks_SYS)
            pks_SYSmax(i)=0;
        else
            [pks_SYSmax(i), index2(i)] = max(pks_SYS);
        end
        
        if isempty(pks_S2)
            pks_S2max(i)=0;
        else
            [pks_S2max(i), index3(i)] = max(pks_S2);
        end
        
        if isempty(pks_DIA)
            pks_DIAmax(i)=0;
            
        else
            
            [pks_DIAmax(i), index4(i)] = max(pks_DIA);
        end
        
    end
    
    if i>1
        % PEAK FREQUENCY
        [pks_S1,~] = findpeaks(2*abs(fft_S1_NFFT(1:NFFT1/2+1)));
        [pks_SYS,~] = findpeaks(2*abs(fft_SYS_NFFT(1:NFFT2/2+1)));
        [pks_S2,~] = findpeaks(2*abs(fft_S2_NFFT(1:NFFT3/2+1)));
        [pks_DIA,~] = findpeaks(2*abs(fft_DIA_NFFT(1:NFFT4/2+1)));
        
        
        if isempty(pks_S1)
            pks_S1max(i)=0;
        else
            [pks_S1max(i), index1(i)] = max(pks_S1);
        end
        
        if isempty(pks_SYS)
            pks_SYSmax(i)=0;
        else
            [pks_SYSmax(i), index2(i)] = max(pks_SYS);
        end
        
        if isempty(pks_S2)
            pks_S2max(i)=0;
        else
            [pks_S2max(i), index3(i)] = max(pks_S2);
        end
        
        if isempty(pks_DIA)
            pks_DIAmax(i)=0;
        else
            [pks_DIAmax(i), index4(i)] = max(pks_DIA);
        end
        
    end
    % BW (SPECTRAL)
%     [~,w1] = periodogram(S1);
%     BWs_S1(i)=bandwidth (w1);
    
    [~,w2] = periodogram(SYS);
    BWs_SYS(i)=bandwidth (w2);
    
%     [~,w3] = periodogram(S2);
%     BWs_S2(i)=bandwidth (w3);
%     
%     [~,w4] = periodogram(DIA);
%     BWs_DIA(i)=bandwidth (w4);
    
    
    % Q-FACTOR (SPECTRAL)
%     Qs1(i)= pks_S1max(i)/BWs_S1(i);
    Qsys(i)= pks_SYSmax(i)/(BWs_SYS(i));
%     Qs2(i)= pks_S2max(i)/(BWs_S2(i));
%     Qdia(i)= pks_DIAmax(i)/(BWs_DIA(i));
    
    % SHANON ENTROPY
    % Compute Shannon entropy of x.
    Shannon_S1(i) = wentropy(S1,'shannon');
    Shannon_SYS(i) = wentropy(SYS,'shannon');
%     Shannon_S2(i) = wentropy(S2,'shannon');
    Shannon_DIA(i) = wentropy(DIA,'shannon');

    % variance
%     var_S1(i)=var(S1);
    var_SYS(i)=var(SYS);
%     var_S2(i)=var(S2);
    var_DIA(i)=var(DIA);
    
    s1I = A(i,1):A(i,2);
    sysI= A(i,2):A(i,3);
    s2I=  A(i,3):A(i,4);
    diaI= A(i,4):A(i+1,1);
    
%     Shannon_A5_s1(i)=wentropy(rA{5}(s1I),'shannon');
%     Shannon_A5_sys(i)=wentropy(rA{5}(sysI),'shannon');
%     Shannon_A5_s2(i)=wentropy(rA{5}(s2I),'shannon');
%     Shannon_A5_dia(i)=wentropy(rA{5}(diaI),'shannon');
%     
%     Shannon_D5_s1(i)=wentropy(rD{5}(s1I),'shannon');
%     Shannon_D5_sys(i)=wentropy(rD{5}(sysI),'shannon');
%     Shannon_D5_s2(i)=wentropy(rD{5}(s2I),'shannon');
%     Shannon_D5_dia(i)=wentropy(rD{5}(diaI),'shannon');
    
%     Shannon_D4_s1(i)=wentropy(rD{4}(s1I),'shannon');
    Shannon_D4_sys(i)=wentropy(rD{4}(sysI),'shannon');
    Shannon_D4_s2(i)=wentropy(rD{4}(s2I),'shannon');
%     Shannon_D4_dia(i)=wentropy(rD{4}(diaI),'shannon');
    
%     Shannon_D3_s1(i)=wentropy(rD{3}(s1I),'shannon');
    Shannon_D3_sys(i)=wentropy(rD{3}(sysI),'shannon');
%     Shannon_D3_s2(i)=wentropy(rD{3}(s2I),'shannon');
%     Shannon_D3_dia(i)=wentropy(rD{3}(diaI),'shannon');
    
    Shannon_D2_s1(i)=wentropy(rD{2}(s1I),'shannon');
    Shannon_D2_sys(i)=wentropy(rD{2}(sysI),'shannon');
    Shannon_D2_s2(i)=wentropy(rD{2}(s2I),'shannon');
    Shannon_D2_dia(i)=wentropy(rD{2}(diaI),'shannon');
    
%     Shannon_D1_s1(i)=wentropy(rD{1}(s1I),'shannon');
%     Shannon_D1_sys(i)=wentropy(rD{1}(sysI),'shannon');
%     Shannon_D1_s2(i)=wentropy(rD{1}(s2I),'shannon');
%     Shannon_D1_dia(i)=wentropy(rD{1}(diaI),'shannon');
end
% mean value of the interval ratios between systole and RR in each heart beat
features.m_Ratio_SysRR   = mean(R_SysRR);
% SD value of the interval ratios between systole and RR in each heart beat
features.sd_Ratio_SysRR  = std(R_SysRR);
% mean value of the interval ratios between diastole and RR in each heart beat
features.m_Ratio_DiaRR   = mean(R_DiaRR);
% SD value of the interval ratios between diastole and RR in each heart beat
features.sd_Ratio_DiaRR  = std(R_DiaRR);
% mean value of the interval ratios between systole and diastole in each heart beat
features.m_Ratio_SysDia  = mean(R_SysDia);
% SD value of the interval ratios between systole and diastole in each heart beat
features.sd_Ratio_SysDia = std(R_SysDia);

% kurtosis
%features.mSK_S1 = mean(SK_S1);
features.sdSK_S1 = std(SK_S1);
% features.mSK_Sys = mean(SK_Sys);
% features.sdSK_Sys = std(SK_Sys);
% features.mSK_S2 = mean(SK_S2);
% features.sdSK_S2 = std(SK_S2);
% features.mSK_Dia = mean(SK_Dia);
% features.sdSK_Dia = std(SK_Dia);

% skewness
features.mKU_S1 = mean(KU_S1);
% features.sdKU_S1 = std(KU_S1);
% features.mKU_Sys = mean(KU_Sys);
% features.sdKU_Sys = std(KU_Sys);
% features.mKU_S2 = mean(KU_S2);
features.sdKU_S2 = std(KU_S2);
% features.mKU_Dia = mean(KU_Dia);
% features.sdKU_Dia = std(KU_Dia);


% avoid the flat line signal
% indx_sys = find(P_SysS1>0 & P_SysS1<100);
% if length(indx_sys)>1
%     % mean value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
%     m_Amp_SysS1  = mean(P_SysS1(indx_sys));
%     % SD value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
%     sd_Amp_SysS1 = std(P_SysS1(indx_sys));
% else
%     m_Amp_SysS1  = 0;
%     sd_Amp_SysS1 = 0;
% end
% indx_dia = find(P_DiaS2>0 & P_DiaS2<100);
% if length(indx_dia)>1
%     % mean value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
%     m_Amp_DiaS2  = mean(P_DiaS2(indx_dia));
%     % SD value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
%     sd_Amp_DiaS2 = std(P_DiaS2(indx_dia));
% else
%     m_Amp_DiaS2  = 0;
%     sd_Amp_DiaS2 = 0;
% end

% Frequency features (potes)
NFFT = 256;
f = (0:NFFT/2-1)/(NFFT/2)*500;
freq_range = [25,45;45,65;65,85;85,105;105,125;125,150;150,200;200,300;300,500];
p_S1  = nan(size(A,1)-1,NFFT/2);
p_Sys = nan(size(A,1)-1,NFFT/2);
p_S2  = nan(size(A,1)-1,NFFT/2);
p_Dia = nan(size(A,1)-1,NFFT/2);
for row=1:size(A,1)-1
    s1 = PCG(A(row,1):A(row,2));
    s1 = s1.*hamming(length(s1));
    Ft = fft(s1,NFFT);
    p_S1(row,:) = abs(Ft(1:NFFT/2));
    
    sys = PCG(A(row,2):A(row,3));
    sys = sys.*hamming(length(sys));
    Ft  = fft(sys,NFFT);
    p_Sys(row,:) = abs(Ft(1:NFFT/2));
    
    s2 = PCG(A(row,3):A(row,4));
    s2 = s2.*hamming(length(s2));
    Ft = fft(s2,NFFT);
    p_S2(row,:) = abs(Ft(1:NFFT/2));
    
    dia = PCG(A(row,4):A(row+1,1));
    dia = dia.*hamming(length(dia));
    Ft  = fft(dia,NFFT);
    p_Dia(row,:) = abs(Ft(1:NFFT/2));
end
P_S1 = nan(1,size(freq_range,1));
P_Sys = nan(1,size(freq_range,1));
P_S2 = nan(1,size(freq_range,1));
P_Dia = nan(1,size(freq_range,1));
for bin=1:size(freq_range,1)
    idx = (f>=freq_range(bin,1)) & (f<freq_range(bin,2));
    P_S1(1,bin) = median(median(p_S1(:,idx)));
    P_Sys(1,bin) = median(median(p_Sys(:,idx)));
    P_S2(1,bin) = median(median(p_S2(:,idx)));
    P_Dia(1,bin) = median(median(p_Dia(:,idx)));
end

features.P_S1=P_S1;
features.P_Sys=P_Sys;
features.P_S2=P_S2;
features.P_Dia=P_Dia;
% features.m_Amp_SysS1=m_Amp_SysS1;  
% features.sd_Amp_SysS1=sd_Amp_SysS1;
% features.m_Amp_DiaS2=m_Amp_DiaS2; 
% features.sd_Amp_DiaS2=sd_Amp_DiaS2;

%% homsi
% % mean value of zero crossing rate for S1, Sys, S2 & Dia
% features.m_Ratio_zcr_S1=mean(zcr_S1);
% features.m_Ratio_zcr_SYS=mean(zcr_SYS);
% features.m_Ratio_zcr_S2=mean(zcr_S2);
% features.m_Ratio_zcr_DIA=mean(zcr_DIA);
% % mean values of duration of S1, Sys, S2 & Dia
% features.m_T_S1=mean(T_S1);
% features.m_T_Sys=mean(T_Sys);
% features.m_T_S2=mean(T_S2);
% features.m_T_Dia=mean(T_Dia);
% mean values of Max
% features.m_Max_S1=mean(Max_S1);
features.m_Max_Sys=mean(Max_Sys);
% features.m_Max_S2=mean(Max_S2);
% features.m_Max_DIA=mean(Max_DIA);
% mean values
% features.m_mean_S1=mean(mean_S1);
% features.m_mean_SYS=mean(mean_SYS);
% features.m_mean_S2=mean(mean_S2);
% features.m_mean_DIA=mean(mean_DIA);
% RMS
%features.m_rms_S1=mean(rms_S1);
features.m_rms_SYS=mean(rms_SYS);
%features.m_rms_S2=mean(rms_S2);
features.m_rms_DIA=mean(rms_DIA);
% mean of Total power in time domain
%features.m_ptd_S1=mean(ptd_S1);
features.m_ptd_SYS=mean(ptd_SYS);
%features.m_ptd_S2=mean(pfd_S2);
features.m_ptd_DIA=mean(ptd_DIA);
% mean of Total power in frequency domain
%features.m_pfd_S1=mean(pfd_S1);
features.m_pfd_SYS=mean(pfd_SYS);
%features.m_pfd_S2=mean(pfd_S2);
features.m_pfd_DIA=mean(pfd_DIA);
% mean Sample Entropy
% features.m_saen_S1=mean(saen_S1);
% features.m_saen_SYS=mean(saen_SYS);
% features.m_saen_S2=mean(saen_S2);
% features.m_saen_DIA=mean(saen_DIA);
% mean values of BW (SPECTRAL)
%features.m_BWs_S1=mean(BWs_S1);
features.m_BWs_SYS=mean(BWs_SYS);
%features.m_BWs_S2=mean(BWs_S2);
%features.m_BWs_DIA=mean(BWs_DIA);
% mean values of Q-FACTOR (SPECTRAL)
%features.m_Qs1=mean(Qs1);
features.m_Qsys=mean(Qsys);
%features.m_Qs2=mean(Qs2);
%features.m_Qdia=mean(Qdia);
% mean values of Shannon entropy
features.m_Shannon_S1=mean(Shannon_S1);
features.m_Shannon_SYS=mean(Shannon_SYS);
%features.m_Shannon_S2=mean(Shannon_S2);
%features.m_Shannon_DIA=mean(Shannon_DIA);
% mean values of variance
%features.m_var_S1=mean(var_S1);
features.m_var_SYS=mean(var_SYS);
%features.m_var_S2=mean(var_S2);
features.m_var_DIA=mean(var_DIA);

% features.m_Shannon_A5_s1 = mean(Shannon_A5_s1);
% features.m_Shannon_A5_sys = mean(Shannon_A5_sys);
% features.m_Shannon_A5_s2 = mean(Shannon_A5_s2);
% features.m_Shannon_A5_dia = mean(Shannon_A5_dia);
% 
% features.m_Shannon_D5_s1 = mean(Shannon_D5_s1);
% features.m_Shannon_D5_sys = mean(Shannon_D5_sys);
% features.m_Shannon_D5_s2 = mean(Shannon_D5_s2);
% features.m_Shannon_D5_dia = mean(Shannon_D5_dia);
% 
% features.m_Shannon_D4_s1 = mean(Shannon_D4_s1);
features.m_Shannon_D4_sys = mean(Shannon_D4_sys);
features.m_Shannon_D4_s2 = mean(Shannon_D4_s2);
%features.m_Shannon_D4_dia = mean(Shannon_D4_dia);

%features.m_Shannon_D3_s1 = mean(Shannon_D3_s1);
features.m_Shannon_D3_sys = mean(Shannon_D3_sys);
%features.m_Shannon_D3_s2 = mean(Shannon_D3_s2);
%features.m_Shannon_D3_dia = mean(Shannon_D3_dia);

% features.m_Shannon_D2_s1 = mean(Shannon_D2_s1);
% features.m_Shannon_D2_sys = mean(Shannon_D2_sys);
% features.m_Shannon_D2_s2 = mean(Shannon_D2_s2);
% features.m_Shannon_D2_dia = mean(Shannon_D2_dia);
% 
% features.m_Shannon_D1_s1 = mean(Shannon_D1_s1);
% features.m_Shannon_D1_sys = mean(Shannon_D1_sys);
% features.m_Shannon_D1_s2 = mean(Shannon_D1_s2);
% features.m_Shannon_D1_dia = mean(Shannon_D1_dia);

% Wavelet Decomposition
% [c,l] = wavedec(PCG,5,'db1');
%Extracting coef of c to be dealt separetely
% A5=c(1:l(1));
% D5=c((l(1)+1):(l(1)+l(2)));
% D4=c((l(1)+l(2)+1):(l(1)+l(2)+l(3)));
% D3=c((l(1)+l(2)+l(3)+1):(l(1)+l(2)+l(3)+l(4)));
% D2=c((l(1)+l(2)+l(3)+l(4)+1):(l(1)+l(2)+l(3)+l(4)+l(5)));
% D1=c((l(1)+l(2)+l(3)+l(4)+l(5)+1):(l(1)+l(2)+l(3)+l(4)+l(5)+l(6)));
% Calculating Shannon entropy for each D + A5
% features.Shannon_A5=wentropy(A5,'shannon');
% features.Shannon_D5=wentropy(D5,'shannon');
% features.Shannon_D4=wentropy(D4,'shannon');
% features.Shannon_D3=wentropy(D3,'shannon');
% features.Shannon_D2=wentropy(D2,'shannon');
% features.Shannon_D1=wentropy(D1,'shannon');

% Calculate entropies of wavelet reconstructed coefficients over each segment type
s1I = [];
s2I = [];
sysI = [];
diaI = [];
for i=1:size(A,1)-1
    s1I =   [s1I A(i,1):A(i,2)];
    sysI = [sysI A(i,2):A(i,3)];
    s2I =   [s2I A(i,3):A(i,4)];
    diaI = [diaI A(i,4):A(i+1,1)];
end

% features.o_Shannon_A5_s1=wentropy(rA{5}(s1I),'shannon');
% features.o_Shannon_A5_sys=wentropy(rA{5}(sysI),'shannon');
% features.o_Shannon_A5_s2=wentropy(rA{5}(s2I),'shannon');
% features.o_Shannon_A5_dia=wentropy(rA{5}(diaI),'shannon');
% 
% features.o_Shannon_D5_s1=wentropy(rD{5}(s1I),'shannon');
% features.o_Shannon_D5_sys=wentropy(rD{5}(sysI),'shannon');
% features.o_Shannon_D5_s2=wentropy(rD{5}(s2I),'shannon');
% features.o_Shannon_D5_dia=wentropy(rD{5}(diaI),'shannon');
% 
% features.o_Shannon_D4_s1=wentropy(rD{4}(s1I),'shannon');
% features.o_Shannon_D4_sys=wentropy(rD{4}(sysI),'shannon');
features.o_Shannon_D4_s2=wentropy(rD{4}(s2I),'shannon');
features.o_Shannon_D4_dia=wentropy(rD{4}(diaI),'shannon');

% features.o_Shannon_D3_s1=wentropy(rD{3}(s1I),'shannon');
% features.o_Shannon_D3_sys=wentropy(rD{3}(sysI),'shannon');
% features.o_Shannon_D3_s2=wentropy(rD{3}(s2I),'shannon');
% features.o_Shannon_D3_dia=wentropy(rD{3}(diaI),'shannon');
% 
% features.o_Shannon_D2_s1=wentropy(rD{2}(s1I),'shannon');
% features.o_Shannon_D2_sys=wentropy(rD{2}(sysI),'shannon');
% features.o_Shannon_D2_s2=wentropy(rD{2}(s2I),'shannon');
% features.o_Shannon_D2_dia=wentropy(rD{2}(diaI),'shannon');
% 
% features.o_Shannon_D1_s1=wentropy(rD{1}(s1I),'shannon');
% features.o_Shannon_D1_sys=wentropy(rD{1}(sysI),'shannon');
% features.o_Shannon_D1_s2=wentropy(rD{1}(s2I),'shannon');
% features.o_Shannon_D1_dia=wentropy(rD{1}(diaI),'shannon');

%detailed_segmentation_features_S1=splitvars(detailed_segmentation_features_S1); 
detailed_segmentation_features_Sys=splitvars(detailed_segmentation_features_Sys); 
detailed_segmentation_features_S2=splitvars(detailed_segmentation_features_S2); 
detailed_segmentation_features_Dia=splitvars(detailed_segmentation_features_Dia); 

%detailed_segmentation_features_S1_mean=array2table(mean(detailed_segmentation_features_S1{:,:})); 
detailed_segmentation_features_Sys_mean=array2table(mean(detailed_segmentation_features_Sys{:,:}));
detailed_segmentation_features_S2_mean=array2table(mean(detailed_segmentation_features_S2{:,:}));
detailed_segmentation_features_Dia_mean=array2table(mean(detailed_segmentation_features_Dia{:,:}));

%detailed_segmentation_features_S1_mean.Properties.VariableNames= strcat(detailed_segmentation_features_S1.Properties.VariableNames,'_S1'); 
detailed_segmentation_features_Sys_mean.Properties.VariableNames= strcat(detailed_segmentation_features_Sys.Properties.VariableNames,'_Sys'); 
detailed_segmentation_features_S2_mean.Properties.VariableNames= strcat(detailed_segmentation_features_S2.Properties.VariableNames,'_S2'); 
detailed_segmentation_features_Dia_mean.Properties.VariableNames= strcat(detailed_segmentation_features_Dia.Properties.VariableNames,'_Dia'); 

features=horzcat(features, detailed_segmentation_features_Sys_mean, detailed_segmentation_features_S2_mean, detailed_segmentation_features_Dia_mean);   
end







