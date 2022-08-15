function features = extractFeaturesFromHsIntervals(assigned_states,PCG, Fs1, heartRate)
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
%             chengyu.liu@emory.edu
%
%  
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

%% Feature calculation
m_RR        = round(mean(diff(A(:,1))));             % mean value of RR intervals
sd_RR       = round(std(diff(A(:,1))));              % standard deviation (SD) value of RR intervals
mean_IntS1  = round(mean(A(:,2)-A(:,1)));            % mean value of S1 intervals
sd_IntS1    = round(std(A(:,2)-A(:,1)));             % SD value of S1 intervals
mean_IntS2  = round(mean(A(:,4)-A(:,3)));            % mean value of S2 intervals
sd_IntS2    = round(std(A(:,4)-A(:,3)));             % SD value of S2 intervals
mean_IntSys = round(mean(A(:,3)-A(:,2)));            % mean value of systole intervals
sd_IntSys   = round(std(A(:,3)-A(:,2)));             % SD value of systole intervals
mean_IntDia = round(mean(A(2:end,1)-A(1:end-1,4)));  % mean value of diastole intervals
sd_IntDia   = round(std(A(2:end,1)-A(1:end-1,4)));   % SD value of diastole intervals

if (size(A,1)-1)>=60
    mysize=60;
else
    mysize=size(A,1)-1;
end;


% PAW 20160803
N_LEVELS = 5;
WAVELET_NAME = 'db4';
DO_PLOT = false;
[C,L] = wavedec(PCG, N_LEVELS, WAVELET_NAME);

if DO_PLOT
    N_ROWS = N_LEVELS + 2;
    figure;
    set(gcf, 'Name', 'Wavelet decomposition');
    subplot(N_ROWS, 1, 1);
    plot(PCG, 'b');
    ylabel('PCG');
    set(gca, 'XTickLabel', []);
    axis tight;
    ax = axis;
    axesH(1) = gca;
end
 
for iLevel = [N_LEVELS:-1:1]
    rD{iLevel} = wrcoef('d', C, L, WAVELET_NAME, iLevel);
    if DO_PLOT
        subplot(N_ROWS, 1, N_LEVELS-iLevel+3);
        plot(rD{iLevel}, 'k');
        axesH(N_LEVELS-iLevel+3) = gca;
        ylabel(sprintf('D%d', iLevel));
        if iLevel > 1
            set(gca, 'XTickLabel', []);
        end
        axis(ax);
    end
    rA{iLevel} = wrcoef('a', C, L, WAVELET_NAME, iLevel);
end

if DO_PLOT
    subplot(N_ROWS, 1, 2);
    plot(rA{N_LEVELS}, 'k');
    axesH(2) = gca;
    ylabel(sprintf('A%d', N_LEVELS));
    set(gca, 'XTickLabel', []);
    axis(ax);
        
    linkaxes(axesH, 'x');
    zoom xon;
end




for i=1:mysize
%for i=1:size(A,1)-1
    R_SysRR(i)  = (A(i,3)-A(i,2))/(A(i+1,1)-A(i,1))*100;
    R_DiaRR(i)  = (A(i+1,1)-A(i,4))/(A(i+1,1)-A(i,1))*100;
    R_SysDia(i) = R_SysRR(i)/R_DiaRR(i)*100;
    
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
    %% ************* New features ***********************************
    %% Signals of S1, Sys, S2 & Dia
    S1=PCG(A(i,1):A(i,2));      L_S1=size(S1);      LS1=L_S1(1,1);
    SYS=PCG(A(i,2):A(i,3));     L_SYS=size(SYS);    LSYS=L_SYS(1,1);
    S2=PCG(A(i,3):A(i,4));      L_S2=size(S2);      LS2=L_S2(1,1);
    DIA=PCG(A(i,4):A(i+1,1));   L_DIA=size(DIA);    LDIA=L_DIA(1,1);
     %% HEART CYCLE 
     HC=PCG (A(i,1):A(i+1,1)-1); 
    
    %% ******* Zero Crossing
    zcr_S1(i) = sum(abs(diff(S1)>0))/length(S1);
    zcr_SYS(i) = sum(abs(diff(SYS)>0))/length(SYS);
    zcr_S2(i) = sum(abs(diff(S2)>0))/length(S2);
    zcr_DIA(i) = sum(abs(diff(DIA)>0))/length(DIA);
    %% ***** Duration of S1, Sys, S2 & Dia 
    T_S1(i)=A(i,2)-A(i,1);
    T_Sys(i)=A(i,3)-A(i,2);
    T_S2(i)=A(i,4)-A(i,3);
    T_Dia(i)=A(i+1,1)-A(i,4);
    %% ***** Maximum 
    Max_S1(i)=max(S1);
    Max_Sys(i)=max(SYS);
    Max_S2(i)=max(S2);
    Max_DIA(i)=max(DIA);
    %%********* Mean
    mean_S1(i)=mean(S1);
    mean_SYS(i)=mean(SYS);
    mean_S2(i)=mean(S2);
    mean_DIA(i)=mean(DIA);
    %%********** RMS
    rms_S1(i)=rms(S1);
    rms_SYS(i)=rms(SYS);
    rms_S2(i)=rms(S2);
    rms_DIA(i)=rms(DIA);
    %% kurtosis
    k_S1(i)=kurtosis(S1);
    k_SYS(i)=kurtosis(SYS);
    k_S2(i)=kurtosis(S2);
    k_DIA(i)=kurtosis(DIA);
    %% Total Power in Time Domain
    ptd_S1(i)=(norm(S1)^2)/LS1;
    ptd_SYS(i)=(norm(SYS)^2)/LSYS;
    ptd_S2(i)=(norm(S2)^2)/LS2;
    ptd_DIA(i)=(norm(DIA)^2)/LDIA;
    %% Total Power in Frequency Domain
    fft_S1=fft(S1);
    fft_SYS=fft(SYS);
    fft_S2=fft(S2);
    fft_DIA=fft(DIA);
    pfd_S1(i)=sum(fft_S1.*conj(fft_S1))/(LS1^2);
    pfd_SYS(i)=sum(fft_SYS.*conj(fft_SYS))/(LSYS^2);
    pfd_S2(i)=sum(fft_S2.*conj(fft_S2))/(LS2^2);
    pfd_DIA(i)=sum(fft_DIA.*conj(fft_DIA))/(LDIA^2);
    
    %% Sample Entropy
    saen_S1(i)=SampEn(2, 0.2*std(S1), S1, 1);
    saen_SYS(i)=SampEn(2, 0.2*std(SYS), SYS, 1);
    saen_S2(i)=SampEn(2, 0.2*std(S2), S2, 1);
    saen_DIA(i)=SampEn(2, 0.2*std(DIA), DIA, 1);
    
    
    
    
    %% FRECUENCY DOMAIN 
    
    NFFT1 = 2^nextpow2(LS1); % Next power of 2 from length of S1
    fft_S1_NFFT = fft(S1,NFFT1)/LS1;
    f1=Fs1/2*linspace(0,1,NFFT1/2+1);
    %plot(f1,2*abs(fft_S1_NFFT(1:NFFT1/2+1)))
    % title('Single-Sided Amplitude Spectrum of S1')
    %xlabel('Frequency (Hz)')
    %ylabel('|S1(f1)|')
     
    NFFT2 = 2^nextpow2(LSYS); % Next power of 2 from length of SYS
    fft_SYS_NFFT = fft(SYS,NFFT2)/LSYS;
    f2=Fs1/2*linspace(0,1,NFFT2/2+1);
    
    NFFT3 = 2^nextpow2(LS2); % Next power of 2 from length of SYS
    fft_S2_NFFT = fft(S2,NFFT3)/LS2;
    f3=Fs1/2*linspace(0,1,NFFT3/2+1);
    
    NFFT4 = 2^nextpow2(LDIA); % Next power of 2 from length of SYS
    fft_DIA_NFFT = fft(DIA,NFFT4)/LDIA;
    f4=Fs1/2*linspace(0,1,NFFT4/2+1);
    
    %% Elaborated by Natasha Medina & Andrea Quintana
    %% conditional
    if i==1
         [pks_S1,locsS1] = findpeaks(2*abs(fft_S1_NFFT(1:NFFT1/2+1)));
         [pks_SYS,locsSYS] = findpeaks(2*abs(fft_SYS_NFFT(1:NFFT2/2+1)));
         [pks_S2,locsS2] = findpeaks(2*abs(fft_S2_NFFT(1:NFFT3/2+1)));
         [pks_DIA,locsDIA] = findpeaks(2*abs(fft_DIA_NFFT(1:NFFT4/2+1)));
       
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
    %% PEAK FREQUENCY 
    [pks_S1,locsS1] = findpeaks(2*abs(fft_S1_NFFT(1:NFFT1/2+1)));
    [pks_SYS,locsSYS] = findpeaks(2*abs(fft_SYS_NFFT(1:NFFT2/2+1)));
    [pks_S2,locsS2] = findpeaks(2*abs(fft_S2_NFFT(1:NFFT3/2+1)));
    [pks_DIA,locsDIA] = findpeaks(2*abs(fft_DIA_NFFT(1:NFFT4/2+1)));
   
    
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
    
    %% Total Harmonic Distortion (THD)
 
     %%THD_S1(i)=thd(2*abs(fft_S1_NFFT(1:NFFT1/2+1)));
     %%THD_SYS(i)=thd(2*abs(fft_SYS_NFFT(1:NFFT2/2+1)));
     %%THD_S2(i)=thd(2*abs(fft_S2_NFFT(1:NFFT3/2+1)));
     %%THD_DIA(i)=thd(2*abs(fft_DIA_NFFT(1:NFFT4/2+1)));
    
    %% BW (SPECTRAL) 
    [pxxS1,w1] = periodogram(S1);
    BWs_S1(i)=bandwidth (w1); 
    
    [pxxSYS,w2] = periodogram(SYS);
    BWs_SYS(i)=bandwidth (w2);
    
    [pxxS2,w3] = periodogram(S2);
    BWs_S2(i)=bandwidth (w3);
    
    [pxxDIA,w4] = periodogram(DIA);
    BWs_DIA(i)=bandwidth (w4);
    
    
    %% Q-FACTOR (SPECTRAL)
    Qs1(i)= pks_S1max(i)/BWs_S1(i);
    Qsys(i)= pks_SYSmax(i)/(BWs_SYS(i));
    Qs2(i)= pks_S2max(i)/(BWs_S2(i));
    Qdia(i)= pks_DIAmax(i)/(BWs_DIA(i));
    
    %% SHANON ENTROPY 
    
    % Compute Shannon entropy of x. 
     Shannon_S1(i) = wentropy(S1,'shannon');
     Shannon_SYS(i) = wentropy(SYS,'shannon');
     Shannon_S2(i) = wentropy(S2,'shannon');
     Shannon_DIA(i) = wentropy(DIA,'shannon');
    
    %% Skewness
    skew_S1(i)=skewness(S1);
    skew_SYS(i)=skewness(SYS);
    skew_S2(i)=skewness(S2);
    skew_DIA(i)=skewness(DIA);
    
    %% variance
    var_S1(i)=var(S1);
    var_SYS(i)=var(SYS);
    var_S2(i)=var(S2);
    var_DIA(i)=var(DIA);
    
    s1I = A(i,1):A(i,2);
    sysI= A(i,2):A(i,3);
    s2I=  A(i,3):A(i,4);
    diaI= A(i,4):A(i+1,1);
    % PAW 20160803
    Shannon_A5_s1(i)=wentropy(rA{5}(s1I),'shannon');
    Shannon_A5_sys(i)=wentropy(rA{5}(sysI),'shannon');
    Shannon_A5_s2(i)=wentropy(rA{5}(s2I),'shannon');
    Shannon_A5_dia(i)=wentropy(rA{5}(diaI),'shannon');
    
    Shannon_D5_s1(i)=wentropy(rD{5}(s1I),'shannon');
    Shannon_D5_sys(i)=wentropy(rD{5}(sysI),'shannon');
    Shannon_D5_s2(i)=wentropy(rD{5}(s2I),'shannon');
    Shannon_D5_dia(i)=wentropy(rD{5}(diaI),'shannon');
    
    Shannon_D4_s1(i)=wentropy(rD{4}(s1I),'shannon');
    Shannon_D4_sys(i)=wentropy(rD{4}(sysI),'shannon');
    Shannon_D4_s2(i)=wentropy(rD{4}(s2I),'shannon');
    Shannon_D4_dia(i)=wentropy(rD{4}(diaI),'shannon');
    
    Shannon_D3_s1(i)=wentropy(rD{3}(s1I),'shannon');
    Shannon_D3_sys(i)=wentropy(rD{3}(sysI),'shannon');
    Shannon_D3_s2(i)=wentropy(rD{3}(s2I),'shannon');
    Shannon_D3_dia(i)=wentropy(rD{3}(diaI),'shannon');
    
    Shannon_D2_s1(i)=wentropy(rD{2}(s1I),'shannon');
    Shannon_D2_sys(i)=wentropy(rD{2}(sysI),'shannon');
    Shannon_D2_s2(i)=wentropy(rD{2}(s2I),'shannon');
    Shannon_D2_dia(i)=wentropy(rD{2}(diaI),'shannon');
    
    Shannon_D1_s1(i)=wentropy(rD{1}(s1I),'shannon');
    Shannon_D1_sys(i)=wentropy(rD{1}(sysI),'shannon');
    Shannon_D1_s2(i)=wentropy(rD{1}(s2I),'shannon');
    Shannon_D1_dia(i)=wentropy(rD{1}(diaI),'shannon');    
end

m_Ratio_SysRR   = mean(R_SysRR);  % mean value of the interval ratios between systole and RR in each heart beat
sd_Ratio_SysRR  = std(R_SysRR);   % SD value of the interval ratios between systole and RR in each heart beat
m_Ratio_DiaRR   = mean(R_DiaRR);  % mean value of the interval ratios between diastole and RR in each heart beat
sd_Ratio_DiaRR  = std(R_DiaRR);   % SD value of the interval ratios between diastole and RR in each heart beat
m_Ratio_SysDia  = mean(R_SysDia); % mean value of the interval ratios between systole and diastole in each heart beat
sd_Ratio_SysDia = std(R_SysDia);  % SD value of the interval ratios between systole and diastole in each heart beat
%% ******************* New Features ****************************
%% mean value of zero crossing rate for S1, Sys, S2 & Dia
m_Ratio_zcr_S1=nanmean(zcr_S1);  
m_Ratio_zcr_SYS=nanmean(zcr_SYS);
m_Ratio_zcr_S2=nanmean(zcr_S2);
m_Ratio_zcr_DIA=nanmean(zcr_DIA);
%% mean values of duration of S1, Sys, S2 & Dia
m_T_S1=nanmean(T_S1);
m_T_Sys=nanmean(T_Sys);
m_T_S2=nanmean(T_S2);
m_T_Dia=nanmean(T_Dia);
%% mean values of Max 
m_Max_S1=nanmean(Max_S1); 
m_Max_Sys=nanmean(Max_Sys);
m_Max_S2=nanmean(Max_S2);
m_Max_DIA=nanmean(Max_DIA);
%% mean values 
m_mean_S1=nanmean(mean_S1); 
m_mean_SYS=nanmean(mean_SYS);
m_mean_S2=nanmean(mean_S2);
m_mean_DIA=nanmean(mean_DIA);
%% RMS
m_rms_S1=nanmean(rms_S1);
m_rms_SYS=nanmean(rms_SYS);
m_rms_S2=nanmean(rms_S2);
m_rms_DIA=nanmean(rms_DIA);
%% kurtosis
m_k_S1=nanmean(k_S1); 
m_k_SYS=nanmean(k_SYS);
m_k_S2=nanmean(k_S2);
m_k_DIA=nanmean(k_DIA);
%% mean of Total power in time domain
m_ptd_S1=nanmean(ptd_S1);
m_ptd_SYS=nanmean(ptd_SYS);
m_ptd_S2=nanmean(pfd_S2);
m_ptd_DIA=nanmean(ptd_DIA);
%% mean of Total power in frequency domain
m_pfd_S1=nanmean(pfd_S1);  
m_pfd_SYS=nanmean(pfd_SYS);
m_pfd_S2=nanmean(pfd_S2);
m_pfd_DIA=nanmean(pfd_DIA);
%% mean Sample Entropy 
m_saen_S1=nanmean(saen_S1);
m_saen_SYS=nanmean(saen_SYS);
m_saen_S2=nanmean(saen_S2);
m_saen_DIA=nanmean(saen_DIA);
%% mean Total Harmonic Distortion (THD)
%m_THD_S1=nanmean (THD_S1);
%m_THD_SYS=nanmean (THD_SYS);
%m_THD_S2=nanmean (THD_S2);
%m_THD_DIA=nanmean (THD_DIA);
 %% mean values of BW (SPECTRAL) 
m_BWs_S1=nanmean(BWs_S1);
m_BWs_SYS=nanmean(BWs_SYS);
m_BWs_S2=nanmean(BWs_S2); 
m_BWs_DIA=nanmean(BWs_DIA);
 %% mean values of Q-FACTOR (SPECTRAL)
 m_Qs1=nanmean(Qs1);
 m_Qsys=nanmean(Qsys);
 m_Qs2=nanmean(Qs2);
 m_Qdia=nanmean(Qdia);
 % mean values of Shannon entropy
 m_Shannon_S1=nanmean(Shannon_S1);
 m_Shannon_SYS=nanmean(Shannon_SYS);
 m_Shannon_S2=nanmean(Shannon_S2);
 m_Shannon_DIA=nanmean(Shannon_DIA);
 % mean values of Skewness
 m_skew_S1=nanmean(skew_S1);
 m_skew_SYS=nanmean(skew_SYS);
 m_skew_S2=nanmean(skew_S2);
 m_skew_DIA=nanmean(skew_DIA);
 % mean values of variance
 m_var_S1=nanmean(var_S1);
 m_var_SYS=nanmean(var_SYS);
 m_var_S2=nanmean(var_S2);
 m_var_DIA=nanmean(var_DIA);
 
indx_sys = find(P_SysS1>0 & P_SysS1<100);   % avoid the flat line signal
if length(indx_sys)>1
    m_Amp_SysS1  = mean(P_SysS1(indx_sys)); % mean value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
    sd_Amp_SysS1 = std(P_SysS1(indx_sys));  % SD value of the mean absolute amplitude ratios between systole period and S1 period in each heart beat
else
    m_Amp_SysS1  = 0;
    sd_Amp_SysS1 = 0;
end
indx_dia = find(P_DiaS2>0 & P_DiaS2<100);
if length(indx_dia)>1
    m_Amp_DiaS2  = mean(P_DiaS2(indx_dia)); % mean value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
    sd_Amp_DiaS2 = std(P_DiaS2(indx_dia));  % SD value of the mean absolute amplitude ratios between diastole period and S2 period in each heart beat
else
    m_Amp_DiaS2  = 0;
    sd_Amp_DiaS2 = 0;
end


m_Shannon_A5_s1 = nanmean(Shannon_A5_s1);
m_Shannon_A5_sys = nanmean(Shannon_A5_sys);
m_Shannon_A5_s2 = nanmean(Shannon_A5_s2);
m_Shannon_A5_dia = nanmean(Shannon_A5_dia);

m_Shannon_D5_s1 = nanmean(Shannon_D5_s1);
m_Shannon_D5_sys = nanmean(Shannon_D5_sys);
m_Shannon_D5_s2 = nanmean(Shannon_D5_s2);
m_Shannon_D5_dia = nanmean(Shannon_D5_dia);

m_Shannon_D4_s1 = nanmean(Shannon_D4_s1);
m_Shannon_D4_sys = nanmean(Shannon_D4_sys);
m_Shannon_D4_s2 = nanmean(Shannon_D4_s2);
m_Shannon_D4_dia = nanmean(Shannon_D4_dia);

m_Shannon_D3_s1 = nanmean(Shannon_D3_s1);
m_Shannon_D3_sys = nanmean(Shannon_D3_sys);
m_Shannon_D3_s2 = nanmean(Shannon_D3_s2);
m_Shannon_D3_dia = nanmean(Shannon_D3_dia);

m_Shannon_D2_s1 = nanmean(Shannon_D2_s1);
m_Shannon_D2_sys = nanmean(Shannon_D2_sys);
m_Shannon_D2_s2 = nanmean(Shannon_D2_s2);
m_Shannon_D2_dia = nanmean(Shannon_D2_dia);

m_Shannon_D1_s1 = nanmean(Shannon_D1_s1);
m_Shannon_D1_sys = nanmean(Shannon_D1_sys);
m_Shannon_D1_s2 = nanmean(Shannon_D1_s2);
m_Shannon_D1_dia = nanmean(Shannon_D1_dia);


if DO_PLOT
    figure;
    set(gcf, 'Name', 'Mean Wavelet entropy');
    compStr = {'A5', 'D5', 'D4', 'D3', 'D2', 'D1'};
    segStr = {'s1', 'sys', 's2', 'dia'};
    barVec = [];
    barVecStr = {};
    for iComp = 1:length(compStr)
        for jSeg = 1:length(segStr)
            barVecStr{end+1} = sprintf('%s_%s', compStr{iComp}, segStr{jSeg});
            barVec(iComp, jSeg) = eval(sprintf('m_Shannon_%s', barVecStr{end}));
        end
    end
    bar(barVec');
    set(gca, 'XTickLabel', segStr);
    legend(compStr);
    ylabel('Shannon entropy');
    
    figure;
    set(gcf, 'Name', 'Wavelet entropy distributions');
    subplot(4, 6, 1); hist(Shannon_A5_s1); xlabel('A5-s1');
    subplot(4, 6, 7); hist(Shannon_A5_sys); xlabel('A5-sys');
    subplot(4, 6, 13); hist(Shannon_A5_s2); xlabel('A5-s2');
    subplot(4, 6, 19); hist(Shannon_A5_dia); xlabel('A5-dia');
    
    subplot(4, 6, 2); hist(Shannon_D5_s1); xlabel('D5-s1');
    subplot(4, 6, 8); hist(Shannon_D5_sys); xlabel('D5-sys');
    subplot(4, 6, 14); hist(Shannon_D5_s2); xlabel('D5-s2');
    subplot(4, 6, 20); hist(Shannon_D5_dia); xlabel('D5-dia');
    
    subplot(4, 6, 3); hist(Shannon_D4_s1); xlabel('D4-s1');
    subplot(4, 6, 9); hist(Shannon_D4_sys); xlabel('D4-sys');
    subplot(4, 6, 15); hist(Shannon_D4_s2); xlabel('D4-s2');
    subplot(4, 6, 21); hist(Shannon_D4_dia); xlabel('D4-dia');
    
    subplot(4, 6, 4); hist(Shannon_D3_s1); xlabel('D3-s1');
    subplot(4, 6, 10); hist(Shannon_D3_sys); xlabel('D3-sys');
    subplot(4, 6, 16); hist(Shannon_D3_s2); xlabel('D3-s2');
    subplot(4, 6, 22); hist(Shannon_D3_dia); xlabel('D3-dia');
    
    subplot(4, 6, 5); hist(Shannon_D2_s2); xlabel('D2-s2');
    subplot(4, 6, 11); hist(Shannon_D2_dia); xlabel('D2-dia');
    subplot(4, 6, 17); hist(Shannon_D2_sys); xlabel('D2-sys');
    subplot(4, 6, 23); hist(Shannon_D2_s1); xlabel('D2-s1');
    
    subplot(4, 6, 6); hist(Shannon_D1_s1); xlabel('D1-s1');
    subplot(4, 6, 12); hist(Shannon_D1_sys); xlabel('D1-sys');
    subplot(4, 6, 18); hist(Shannon_D1_s2); xlabel('D1-s2');
    subplot(4, 6, 24); hist(Shannon_D1_dia); xlabel('D1-dia');
end

%% Elaborated by Miguel Hernandez
%PCG1=PCG([1:ceil(size(PCG,1)/2)]);
%% Wavelet Decomposition 
[c,l] = wavedec(PCG,5,'db1');
[Ea Ed]=wenergy(c,l); %%Ea= Energy of A5. Ed= Energy of coef D
%Extracting coef of c to be dealt separetely
A5=c(1:l(1));
D5=c((l(1)+1):(l(1)+l(2)));
D4=c((l(1)+l(2)+1):(l(1)+l(2)+l(3)));
D3=c((l(1)+l(2)+l(3)+1):(l(1)+l(2)+l(3)+l(4)));
D2=c((l(1)+l(2)+l(3)+l(4)+1):(l(1)+l(2)+l(3)+l(4)+l(5)));
D1=c((l(1)+l(2)+l(3)+l(4)+l(5)+1):(l(1)+l(2)+l(3)+l(4)+l(5)+l(6)));
% Calculating Shannon entropy for each D + A5
Shannon_A5=wentropy(A5,'shannon');
Shannon_D5=wentropy(D5,'shannon');
Shannon_D4=wentropy(D4,'shannon');
Shannon_D3=wentropy(D3,'shannon');
Shannon_D2=wentropy(D2,'shannon');
Shannon_D1=wentropy(D1,'shannon');

%
% PAW: added 20160730
%      -calculate entropies of wavelet reconstructed coefficients
%       over each segment type
%

s1I = [];
s2I = [];
sysI = [];
diaI = [];
for i=1:mysize
    s1I =   [s1I A(i,1):A(i,2)];
    sysI = [sysI A(i,2):A(i,3)];
    s2I =   [s2I A(i,3):A(i,4)];
    diaI = [diaI A(i,4):A(i+1,1)];
end

o_Shannon_A5_s1=wentropy(rA{5}(s1I),'shannon');
o_Shannon_A5_sys=wentropy(rA{5}(sysI),'shannon');
o_Shannon_A5_s2=wentropy(rA{5}(s2I),'shannon');
o_Shannon_A5_dia=wentropy(rA{5}(diaI),'shannon');

o_Shannon_D5_s1=wentropy(rD{5}(s1I),'shannon');
o_Shannon_D5_sys=wentropy(rD{5}(sysI),'shannon');
o_Shannon_D5_s2=wentropy(rD{5}(s2I),'shannon');
o_Shannon_D5_dia=wentropy(rD{5}(diaI),'shannon');

o_Shannon_D4_s1=wentropy(rD{4}(s1I),'shannon');
o_Shannon_D4_sys=wentropy(rD{4}(sysI),'shannon');
o_Shannon_D4_s2=wentropy(rD{4}(s2I),'shannon');
o_Shannon_D4_dia=wentropy(rD{4}(diaI),'shannon');

o_Shannon_D3_s1=wentropy(rD{3}(s1I),'shannon');
o_Shannon_D3_sys=wentropy(rD{3}(sysI),'shannon');
o_Shannon_D3_s2=wentropy(rD{3}(s2I),'shannon');
o_Shannon_D3_dia=wentropy(rD{3}(diaI),'shannon');

o_Shannon_D2_s1=wentropy(rD{2}(s1I),'shannon');
o_Shannon_D2_sys=wentropy(rD{2}(sysI),'shannon');
o_Shannon_D2_s2=wentropy(rD{2}(s2I),'shannon');
o_Shannon_D2_dia=wentropy(rD{2}(diaI),'shannon');

o_Shannon_D1_s1=wentropy(rD{1}(s1I),'shannon');
o_Shannon_D1_sys=wentropy(rD{1}(sysI),'shannon');
o_Shannon_D1_s2=wentropy(rD{1}(s2I),'shannon');
o_Shannon_D1_dia=wentropy(rD{1}(diaI),'shannon');

if DO_PLOT
    figure;
    set(gcf, 'Name', 'Overall Wavelet entropy');
    compStr = {'A5', 'D5', 'D4', 'D3', 'D2', 'D1'};
    segStr = {'s1', 'sys', 's2', 'dia'};
    barVec = [];
    barVecStr = {};
    for iComp = 1:length(compStr)
        for jSeg = 1:length(segStr)
            barVecStr{end+1} = sprintf('%s_%s', compStr{iComp}, segStr{jSeg});
            barVec(iComp, jSeg) = eval(sprintf('o_Shannon_%s', barVecStr{end}));
        end
    end
    bar(barVec');
    set(gca, 'XTickLabel', segStr);
    legend(compStr);
    ylabel('Shannon entropy');
end
%
% PAW: end
%

MyMeans = [875.39 43.481 131.188 14.735 105.391 ...
    13.483 195.65 15.789 443.293 32.769...
    22.363 1.71 49.025 2.18 48.308...
    5.078 38.487 10.444 47.723 11.563...
    0.494 0.496 0.495 0.495 131.362...
    195.204 105.277 443.66 0.331 0.103...
    0.217 0.12 -0.02 -0.02 -0.021...
    -0.021 0.154 0.064 0.107 0.06...
    3.98 4.881 4.452 7.007 0.037...
    0.011 0.02 0.009 0.037 0.011...
    0.02 0.009 0.452 0.694 0.481...
    0.68 128 144.153 128 325.668...
    0.001 0 0.001 0 8.839...
    5.442 4.753 10.952 -0.03 -0.061...
    -0.057 -0.022 0.032 0.006 0.015...
    0.005 73.565 -68.39 15.45 38.001...
    56.617 45.63 25.933 3.589 4.37...
    2.168 8.953 3.855 1.134 1.413...
    2.136 3.678 0.609 1.692 1.172...
    1.642 0.364 0.978 0.565 0.419...
    0.17 0.335 0.26 0.092 0.057...
    0.071 0.1 47.285 47.63 26.475...
    96.511 75.211 17.045 24.927 32.175...
    84.156 11.586 37.54 22.951 42.993...
    8.585 24.252 12.922 12.475 4.489...
    9.121 6.526 3.22 1.607 2.26...
    2.637];
features = [m_RR sd_RR  mean_IntS1 sd_IntS1  mean_IntS2 sd_IntS2  mean_IntSys sd_IntSys  mean_IntDia sd_IntDia ...
    m_Ratio_SysRR sd_Ratio_SysRR m_Ratio_DiaRR sd_Ratio_DiaRR m_Ratio_SysDia sd_Ratio_SysDia ...
    m_Amp_SysS1 sd_Amp_SysS1 m_Amp_DiaS2 sd_Amp_DiaS2 ...
    m_Ratio_zcr_S1 m_Ratio_zcr_SYS m_Ratio_zcr_S2 m_Ratio_zcr_DIA ...
    m_T_S1 m_T_Sys m_T_S2 m_T_Dia ...
    m_Max_S1 m_Max_Sys m_Max_S2 m_Max_DIA ...
    m_mean_S1 m_mean_SYS m_mean_S2 m_mean_DIA ...
    m_rms_S1 m_rms_SYS m_rms_S2 m_rms_DIA ...
    m_k_S1 m_k_SYS m_k_S2 m_k_DIA ...
    m_ptd_S1 m_ptd_SYS m_ptd_S2 m_ptd_DIA ...
    m_pfd_S1 m_pfd_SYS m_pfd_S2 m_pfd_DIA ...
    m_saen_S1 m_saen_SYS m_saen_S2 m_saen_DIA ...
    m_BWs_S1 m_BWs_SYS m_BWs_S2 m_BWs_DIA ...
    m_Qs1 m_Qsys m_Qs2 m_Qdia ...
    m_Shannon_S1 m_Shannon_SYS m_Shannon_S2 m_Shannon_DIA ...
    m_skew_S1 m_skew_SYS m_skew_S2 m_skew_DIA ...
    m_var_S1 m_var_SYS m_var_S2 m_var_DIA ...
    heartRate...
    Shannon_A5 Shannon_D5 Shannon_D4 Shannon_D3 Shannon_D2 Shannon_D1 ...
...
    m_Shannon_A5_s1 m_Shannon_A5_sys m_Shannon_A5_s2 m_Shannon_A5_dia ...
    m_Shannon_D5_s1 m_Shannon_D5_sys m_Shannon_D5_s2 m_Shannon_D5_dia ...
    m_Shannon_D4_s1 m_Shannon_D4_sys m_Shannon_D4_s2 m_Shannon_D4_dia ...
    m_Shannon_D3_s1 m_Shannon_D3_sys m_Shannon_D3_s2 m_Shannon_D3_dia ...
    m_Shannon_D2_s1 m_Shannon_D2_sys m_Shannon_D2_s2 m_Shannon_D2_dia ...
    m_Shannon_D1_s1 m_Shannon_D1_sys m_Shannon_D1_s2 m_Shannon_D1_dia ...
...
    o_Shannon_A5_s1 o_Shannon_A5_sys o_Shannon_A5_s2 o_Shannon_A5_dia ...
    o_Shannon_D5_s1 o_Shannon_D5_sys o_Shannon_D5_s2 o_Shannon_D5_dia ...
    o_Shannon_D4_s1 o_Shannon_D4_sys o_Shannon_D4_s2 o_Shannon_D4_dia ...
    o_Shannon_D3_s1 o_Shannon_D3_sys o_Shannon_D3_s2 o_Shannon_D3_dia ...
    o_Shannon_D2_s1 o_Shannon_D2_sys o_Shannon_D2_s2 o_Shannon_D2_dia ...
    o_Shannon_D1_s1 o_Shannon_D1_sys o_Shannon_D1_s2 o_Shannon_D1_dia];


TF = isinf(features);
Zeroos=[features==0];
for i = 1:131
    if TF(i)==1 || Zeroos(i)==1
     features(i)=MyMeans(i);
    end;
end;


if DO_PLOT
    close all;
end



