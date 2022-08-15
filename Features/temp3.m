% Plesinger
%% Load data and resample data
[PCG, Fs1]

ultraFEnvelope = hobalka(PCG,Fs1,400,800);

Fs = 1000;
PCG_resampled      = resample(PCG,Fs,Fs1); 

ultraFEnvelope = resample(ultraFEnvelope,Fs,Fs1);

lfMin = 15;
lfMax = 90;

mfMin = 55;
mfMax = 150;

hfMin = 100;
hfMax = 250;

sfMin = 200;
sfMax = 450;

invalid_window=1000;
invalid_step = 500;
invalid = zeros(length(PCG_resampled),1);

lowFEnvelope = hobalka(PCG_resampled,Fs,lfMin,lfMax);
midFEnvelope = hobalka(PCG_resampled,Fs,mfMin,mfMax);
highFEnvelope = hobalka(PCG_resampled,Fs,hfMin,hfMax);
superFEnvelope = hobalka(PCG_resampled,Fs,sfMin,sfMax);


for s=1:invalid_step:length(PCG_resampled)-invalid_window
    testedBlock = PCG_resampled(s:s+invalid_window);
    hst = hist(testedBlock,10);
    
    hfRatio = sum(highFEnvelope(s:s+invalid_window))/sum(lowFEnvelope(s:s+invalid_window));
    
    invalidityResult =(sum(hst(1:2))+sum(hst(9:10)))>sum(hst(3:8));
    invalidityResult = invalidityResult | sum(hst(1))>mean(hst(2:9)) | sum(hst(10))>mean(hst(2:9));
    invalidityResult = invalidityResult | (max(testedBlock)-min(testedBlock)==0);
    invalidityResult = invalidityResult | (hfRatio>1);
    
    invalid(s:s+invalid_window) = invalidityResult;
end


invalidPart = mean(invalid);
globalHFratio = sum(highFEnvelope)/sum(lowFEnvelope);


%frequency ratio
env_1_40 = hobalka(PCG_resampled,Fs,1,40);
env_60_200 = hobalka(PCG_resampled,Fs,60,200);
frequencyRatio = sum(env_60_200)/sum(env_1_40);

[s1positions, s2positions, extrema, locExtrema, minHistoVal, maxHistoVal,...
    histogram, avgLShape, avgMShape, avgHShape, avgSShape, avgUShape,...
    mdetcorrs, fftcorrs] = getS1S2(PCG_resampled,Fs,lowFEnvelope,midFEnvelope,highFEnvelope,superFEnvelope, ultraFEnvelope, invalid );

rrs = diff(locExtrema(s1positions));

% for statistics
s1s2Widths = locExtrema(s2positions)-locExtrema(s1positions); 

% Extract features
ri=0;

esum = 0.6;
SK2M = 1.7;

for index=1:length(s1positions)
   
   index2 = fix(index);
    
   s1pos = locExtrema(s1positions(index2));
   s2pos = locExtrema(s2positions(index2));
   sysW = s2pos-s1pos; 
   
   if mod(sysW,2)>0
       sysW=sysW+1;
   end
   
   soundW = sysW/4;
   
   startA=round(s1pos-sysW/2);
   startB=round(s1pos-soundW/2);
   startC=round(s1pos+soundW/2);
   startD=round(s1pos+sysW/2);
   startE=round(s2pos-soundW/2);
   startF=round(s2pos+soundW/2);
   endF=round(s2pos+sysW/2);
   endG=round(s2pos+sysW);
   
   if startA<=0
       continue;
   end
   
   if startF>length(PCG_resampled) || endF>length(PCG_resampled)
       break;
   end
   
   if endG>length(PCG_resampled)
       endG = length(PCG_resampled);
   end
   
   ri = ri+1;
   
   startSymS1=s1pos-sysW/2;
   endSymS1=s1pos+sysW/2;
   startSymS2=s2pos-sysW/2;
   endSymS2 = s2pos+sysW/2;
   
   leftS1=fliplr(lowFEnvelope(startSymS1:s1pos)');
   rightS1=lowFEnvelope(s1pos:endSymS1)';
   
   leftS2=fliplr(lowFEnvelope(startSymS2:s2pos)');
   rightS2=lowFEnvelope(s2pos:endSymS2)';
   
   diffS1=abs(leftS1-rightS1);
   diffS2=abs(leftS2-rightS2);
   
   a1 = corrcoef(leftS1,rightS1);
   a2 = corrcoef(leftS2,rightS2);
    
   symmetryS1sum(ri)=sum(diffS1);
   symmetryS1std(ri)=std(diffS1);
   symmetryS1Cor(ri)=a1(1,2);
   
   symmetryS2sum(ri)=sum(diffS2);
   symmetryS2std(ri)=std(diffS2);
   symmetryS2Cor(ri)=a2(1,2);
   
   leftS1=fliplr(midFEnvelope(startSymS1:s1pos)'); 
   rightS1=midFEnvelope(s1pos:endSymS1)';
   leftS2=fliplr(midFEnvelope(startSymS2:s2pos)');
   rightS2=midFEnvelope(s2pos:endSymS2)';
   a1 = corrcoef(leftS1,rightS1);
   a2 = corrcoef(leftS2,rightS2);
   symmetryS1CorMF(ri)=a1(1,2);
   symmetryS2CorMF(ri)=a2(1,2);
   
   leftS1=fliplr(highFEnvelope(startSymS1:s1pos)'); 
   rightS1=highFEnvelope(s1pos:endSymS1)';
   leftS2=fliplr(highFEnvelope(startSymS2:s2pos)');
   rightS2=highFEnvelope(s2pos:endSymS2)';
   a1 = corrcoef(leftS1,rightS1);
   a2 = corrcoef(leftS2,rightS2);
   symmetryS1CorHF(ri)=a1(1,2);
   symmetryS2CorHF(ri)=a2(1,2);
   
   leftS1=fliplr(superFEnvelope(startSymS1:s1pos)'); 
   rightS1=superFEnvelope(s1pos:endSymS1)';
   leftS2=fliplr(superFEnvelope(startSymS2:s2pos)');
   rightS2=superFEnvelope(s2pos:endSymS2)';
   a1 = corrcoef(leftS1,rightS1);
   a2 = corrcoef(leftS2,rightS2);
   symmetryS1CorSF(ri)=a1(1,2);
   symmetryS2CorSF(ri)=a2(1,2);
   
   leftS1=fliplr(ultraFEnvelope(startSymS1:s1pos)'); 
   rightS1=ultraFEnvelope(s1pos:endSymS1)';
   leftS2=fliplr(ultraFEnvelope(startSymS2:s2pos)');
   rightS2=ultraFEnvelope(s2pos:endSymS2)';
   a1 = corrcoef(leftS1,rightS1);
   a2 = corrcoef(leftS2,rightS2);
   symmetryS1CorUF(ri)=a1(1,2);
   symmetryS2CorUF(ri)=a2(1,2);
   
   skewnessRAW(ri)=skewness(PCG_resampled(startSymS1:endSymS2));
   skewnessLOW(ri)=skewness(lowFEnvelope(startSymS1:endSymS2));
   skewnessMID(ri)=skewness(midFEnvelope(startSymS1:endSymS2));
   skewnessHIGH(ri)=skewness(highFEnvelope(startSymS1:endSymS2));
   skewnessSUPER(ri)=skewness(superFEnvelope(startSymS1:endSymS2));
   skewnessULTRA(ri)=skewness(ultraFEnvelope(startSymS1:endSymS2));
   
   skewS1LOW(ri) = skewness(lowFEnvelope(startSymS1:endSymS1));
   skewS1MID(ri) = skewness(midFEnvelope(startSymS1:endSymS1));
   skewS1HIGH(ri) = skewness(highFEnvelope(startSymS1:endSymS1));
   skewS1SUPER(ri) = skewness(superFEnvelope(startSymS1:endSymS1));
   skewS1ULTRA(ri) = skewness(ultraFEnvelope(startSymS1:endSymS1));
   
   skewS2LOW(ri) = skewness(lowFEnvelope(startSymS2:endSymS2));
   skewS2MID(ri) = skewness(midFEnvelope(startSymS2:endSymS2));
   skewS2HIGH(ri) = skewness(highFEnvelope(startSymS2:endSymS2));
   skewS2SUPER(ri) = skewness(superFEnvelope(startSymS2:endSymS2));
   skewS2ULTRA(ri) = skewness(ultraFEnvelope(startSymS2:endSymS2));
   
   dataAsum(ri)=sum(env_60_200(startA:startB));
   dataBsum(ri)=sum(env_60_200(startB:startC));
   dataCsum(ri)=sum(env_60_200(startC:startD));
   dataDsum(ri)=sum(env_60_200(startD:startE));
   dataEsum(ri)=sum(env_60_200(startE:startF));
   dataFsum(ri)=sum(env_60_200(startF:endF));
   dataGsum(ri)=sum(env_60_200(endF:endG));
   
   dataAstd(ri)=std(env_60_200(startA:startB));
   dataBstd(ri)=std(env_60_200(startB:startC));
   dataCstd(ri)=std(env_60_200(startC:startD));
   dataDstd(ri)=std(env_60_200(startD:startE));
   dataEstd(ri)=std(env_60_200(startE:startF));
   dataFstd(ri)=std(env_60_200(startF:endF));
   dataGstd(ri)=std(env_60_200(endF:endG));
   
   dataAmed(ri)=median(env_60_200(startA:startB));
   dataBmed(ri)=median(env_60_200(startB:startC));
   dataCmed(ri)=median(env_60_200(startC:startD));
   dataDmed(ri)=median(env_60_200(startD:startE));
   dataEmed(ri)=median(env_60_200(startE:startF));
   dataFmed(ri)=median(env_60_200(startF:endF));
   dataGmed(ri)=median(env_60_200(endF:endG));

   SdataAsum(ri)=sum(superFEnvelope(startA:startB));
   SdataBsum(ri)=sum(superFEnvelope(startB:startC));
   SdataCsum(ri)=sum(superFEnvelope(startC:startD));
   SdataDsum(ri)=sum(superFEnvelope(startD:startE));
   SdataEsum(ri)=sum(superFEnvelope(startE:startF));
   SdataFsum(ri)=sum(superFEnvelope(startF:endF));
   SdataGsum(ri)=sum(superFEnvelope(endF:endG));
   
   SdataAstd(ri)=std(superFEnvelope(startA:startB));
   SdataBstd(ri)=std(superFEnvelope(startB:startC));
   SdataCstd(ri)=std(superFEnvelope(startC:startD));
   SdataDstd(ri)=std(superFEnvelope(startD:startE));
   SdataEstd(ri)=std(superFEnvelope(startE:startF));
   SdataFstd(ri)=std(superFEnvelope(startF:endF));
   SdataGstd(ri)=std(superFEnvelope(endF:endG));
   
   SdataAmed(ri)=median(superFEnvelope(startA:startB));
   SdataBmed(ri)=median(superFEnvelope(startB:startC));
   SdataCmed(ri)=median(superFEnvelope(startC:startD));
   SdataDmed(ri)=median(superFEnvelope(startD:startE));
   SdataEmed(ri)=median(superFEnvelope(startE:startF));
   SdataFmed(ri)=median(superFEnvelope(startF:endF));
   SdataGmed(ri)=median(superFEnvelope(endF:endG));
   
   UdataAsum(ri)=sum(ultraFEnvelope(startA:startB));
   UdataBsum(ri)=sum(ultraFEnvelope(startB:startC));
   UdataCsum(ri)=sum(ultraFEnvelope(startC:startD));
   UdataDsum(ri)=sum(ultraFEnvelope(startD:startE));
   UdataEsum(ri)=sum(ultraFEnvelope(startE:startF));
   UdataFsum(ri)=sum(ultraFEnvelope(startF:endF));
   UdataGsum(ri)=sum(ultraFEnvelope(endF:endG));
   
   UdataAstd(ri)=std(ultraFEnvelope(startA:startB));
   UdataBstd(ri)=std(ultraFEnvelope(startB:startC));
   UdataCstd(ri)=std(ultraFEnvelope(startC:startD));
   UdataDstd(ri)=std(ultraFEnvelope(startD:startE));
   UdataEstd(ri)=std(ultraFEnvelope(startE:startF));
   UdataFstd(ri)=std(ultraFEnvelope(startF:endF));
   UdataGstd(ri)=std(ultraFEnvelope(endF:endG));
   
   UdataAmed(ri)=median(ultraFEnvelope(startA:startB));
   UdataBmed(ri)=median(ultraFEnvelope(startB:startC));
   UdataCmed(ri)=median(ultraFEnvelope(startC:startD));
   UdataDmed(ri)=median(ultraFEnvelope(startD:startE));
   UdataEmed(ri)=median(ultraFEnvelope(startE:startF));
   UdataFmed(ri)=median(ultraFEnvelope(startF:endF));
   UdataGmed(ri)=median(ultraFEnvelope(endF:endG));
   
   if (startF<length(mdetcorrs))
   jhPosesGlobal = startB:startF;
   jhPosesLocal = jhPosesGlobal-startD;
   
   arg1 = sum(mdetcorrs(jhPosesGlobal)'.*jhPosesLocal.*jhPosesLocal);
   arg2 = sum(mdetcorrs(jhPosesGlobal));
   
   jhParam(ri)=arg1/arg2;
   end
end

Corr_FFT_mean = mean(mdetcorrs);
Corr_FFT_std = std(mdetcorrs);
Corr_FFT_skw = skewness(mdetcorrs);
Corr_FFT_kurt = kurtosis(mdetcorrs);
Corr_FFT_median = median(mdetcorrs);

Corr_FFT_mean_i_std = Corr_FFT_mean/Corr_FFT_std;

jhParam_mean = mean(jhParam);
jhParam_std = std(jhParam);
jhParam_skew = skewness(jhParam);
jhParam_kurt = kurtosis(jhParam);
    
AM_mean = mean(dataAmed);
BM_mean = mean(dataBmed);
CM_mean = mean(dataCmed);
DM_mean = mean(dataDmed);
EM_mean = mean(dataEmed);
FM_mean = mean(dataFmed);
GM_mean = mean(dataGmed);

AD_mean = mean(dataAstd);
BD_mean = mean(dataBstd);
CD_mean = mean(dataCstd);
DD_mean = mean(dataDstd);
ED_mean = mean(dataEstd);
FD_mean = mean(dataFstd);
GD_mean = mean(dataGstd);

AS_mean = mean(dataAsum);
BS_mean = mean(dataBsum);
CS_mean = mean(dataCsum);
DS_mean = mean(dataDsum);
ES_mean = mean(dataEsum);
FS_mean = mean(dataFsum);
GS_mean = mean(dataGsum);

ADiBD = AD_mean/BD_mean;
ADiCD = AD_mean/CD_mean;
CDiDD = CD_mean/DD_mean;
DDiFD = DD_mean/FD_mean;
BDiED = BD_mean/ED_mean;
EDiFD = ED_mean/FD_mean;
DDiED = DD_mean/ED_mean;
ADiGD = AD_mean/GD_mean;

AMiBM = AM_mean/BM_mean;
AMiCM = AM_mean/CM_mean;
BMiEM = BM_mean/EM_mean;
BMiCM = BM_mean/CM_mean;
DMiFM = DM_mean/FM_mean;
DMiEM = DM_mean/EM_mean;
AMiGM = AM_mean/GM_mean;

ACDFMiBEM=(AM_mean+CM_mean+DM_mean+FM_mean)/(BM_mean+EM_mean);

%superF
SAM_mean = mean(SdataAmed);
SBM_mean = mean(SdataBmed);
SCM_mean = mean(SdataCmed);
SDM_mean = mean(SdataDmed);
SEM_mean = mean(SdataEmed);
SFM_mean = mean(SdataFmed);
SGM_mean = mean(SdataGmed);

SAD_mean = mean(SdataAstd);
SBD_mean = mean(SdataBstd);
SCD_mean = mean(SdataCstd);
SDD_mean = mean(SdataDstd);
SED_mean = mean(SdataEstd);
SFD_mean = mean(SdataFstd);
SGD_mean = mean(SdataGstd);

SAS_mean = mean(SdataAsum);
SBS_mean = mean(SdataBsum);
SCS_mean = mean(SdataCsum);
SDS_mean = mean(SdataDsum);
SES_mean = mean(SdataEsum);
SFS_mean = mean(SdataFsum);
SGS_mean = mean(SdataGsum);

SADiBD = SAD_mean/SBD_mean;
SADiCD = SAD_mean/SCD_mean;
SCDiDD = SCD_mean/SDD_mean;
SDDiFD = SDD_mean/SFD_mean;
SBDiED = SBD_mean/SED_mean;
SEDiFD = SED_mean/SFD_mean;
SDDiED = SDD_mean/SED_mean;
SADiGD = SAD_mean/SGD_mean;

SAMiBM = SAM_mean/SBM_mean;
SAMiCM = SAM_mean/SCM_mean;
SBMiEM = SBM_mean/SEM_mean;
SBMiCM = SBM_mean/SCM_mean;
SDMiFM = SDM_mean/SFM_mean;
SDMiEM = SDM_mean/SEM_mean;
SAMiGM = SAM_mean/SGM_mean;

SACDFMiBEM=(SAM_mean+SCM_mean+SDM_mean+SFM_mean)/(SBM_mean+SEM_mean);

%ultraF
UAM_mean = mean(UdataAmed);
UBM_mean = mean(UdataBmed);
UCM_mean = mean(UdataCmed);
UDM_mean = mean(UdataDmed);
UEM_mean = mean(UdataEmed);
UFM_mean = mean(UdataFmed);
UGM_mean = mean(UdataGmed);

UAD_mean = mean(UdataAstd);
UBD_mean = mean(UdataBstd);
UCD_mean = mean(UdataCstd);
UDD_mean = mean(UdataDstd);
UED_mean = mean(UdataEstd);
UFD_mean = mean(UdataFstd);
UGD_mean = mean(UdataGstd);

UAS_mean = mean(UdataAsum);
UBS_mean = mean(UdataBsum);
UCS_mean = mean(UdataCsum);
UDS_mean = mean(UdataDsum);
UES_mean = mean(UdataEsum);
UFS_mean = mean(UdataFsum);
UGS_mean = mean(UdataGsum);

UADiBD = UAD_mean/UBD_mean;
UADiCD = UAD_mean/UCD_mean;
UCDiDD = UCD_mean/UDD_mean;
UDDiFD = UDD_mean/UFD_mean;
UBDiED = UBD_mean/UED_mean;
UEDiFD = UED_mean/UFD_mean;
UDDiED = UDD_mean/UED_mean;
UADiGD = UAD_mean/UGD_mean;

UAMiBM = UAM_mean/UBM_mean;
UAMiCM = UAM_mean/UCM_mean;
UBMiEM = UBM_mean/UEM_mean;
UBMiCM = UBM_mean/UCM_mean;
UDMiFM = UDM_mean/UFM_mean;
UDMiEM = UDM_mean/UEM_mean;
UAMiGM = UAM_mean/UGM_mean;

UACDFMiBEM=(UAM_mean+UCM_mean+UDM_mean+UFM_mean)/(UBM_mean+UEM_mean);
UAMiAM = UAM_mean/AM_mean;
UBMiBM = UBM_mean/BM_mean;
UCMiCM = UCM_mean/CM_mean;
UDMiDM = UDM_mean/DM_mean;
UEMiEM = UEM_mean/EM_mean;
UFMiFM = UFM_mean/FM_mean;
UGMiGM = UGM_mean/GM_mean;

% score calculation
SystoleWidth_std = std(s1s2Widths);
SystoleWidth_mean = mean(s1s2Widths);

RR_std=std(rrs);
RR_count=length(rrs);
RR_mean = mean(rrs);

Diastoles = [];
for d=2:length(s2positions)
    Diastoles(d-1) = locExtrema(s1positions(d))-locExtrema(s2positions(d-1));
end

Diastole_mean = mean(Diastoles);
Diastole_std = std(Diastoles);

upSystoles = s1s2Widths';
upSystoles(end)=[];

SysDiaRatio = Diastoles./upSystoles;

Sys_i_dia_mean = mean(SysDiaRatio);
Sys_i_dia_std = std(SysDiaRatio);

Sym_S1_Corr_mean = mean(symmetryS1Cor);
Sym_S2_Corr_mean = mean(symmetryS2Cor);
Sym_S1_Corr_std = std(symmetryS1Cor);
Sym_S2_Corr_std = std(symmetryS2Cor);
Sym_S1_Sum_mean = mean(symmetryS1sum);
Sym_S2_Sum_mean = mean(symmetryS2sum);
Sym_S1_Sum_std=std(symmetryS1sum);
Sym_S2_Sum_std=std(symmetryS2sum);

Sym_S1_Sum_mean_i_std = Sym_S1_Sum_mean/Sym_S1_Sum_std;
Sym_S2_Sum_mean_i_std = Sym_S2_Sum_mean/Sym_S2_Sum_std;

Sym_S1_Corr_mean_i_std = Sym_S1_Corr_mean/Sym_S1_Corr_std;
Sym_S2_Corr_mean_i_std = Sym_S2_Corr_mean/Sym_S2_Corr_std;

Sym_S1_Corr_MF_mean = mean(symmetryS1CorMF);
Sym_S2_Corr_MF_mean = mean(symmetryS2CorMF);
Sym_S1_Corr_HF_mean = mean(symmetryS1CorHF);
Sym_S2_Corr_HF_mean = mean(symmetryS2CorHF);
Sym_S1_Corr_SF_mean = mean(symmetryS1CorSF);
Sym_S2_Corr_SF_mean = mean(symmetryS2CorSF);
Sym_S1_Corr_UF_mean = mean(symmetryS1CorUF);
Sym_S2_Corr_UF_mean = mean(symmetryS2CorUF);

Sym_S1_Corr_MF_std = std(symmetryS1CorMF);
Sym_S2_Corr_MF_std = std(symmetryS2CorMF);
Sym_S1_Corr_HF_std = std(symmetryS1CorHF);
Sym_S2_Corr_HF_std = std(symmetryS2CorHF);
Sym_S1_Corr_SF_std = std(symmetryS1CorSF);
Sym_S2_Corr_SF_std = std(symmetryS2CorSF);
Sym_S1_Corr_UF_std = std(symmetryS1CorUF);
Sym_S2_Corr_UF_std = std(symmetryS2CorUF);

Sym_S1_Corr_MF_mean_i_std = Sym_S1_Corr_MF_mean/Sym_S1_Corr_MF_std;
Sym_S2_Corr_MF_mean_i_std = Sym_S2_Corr_MF_mean/Sym_S2_Corr_MF_std;
Sym_S1_Corr_HF_mean_i_std = Sym_S1_Corr_HF_mean/Sym_S1_Corr_HF_std;
Sym_S2_Corr_HF_mean_i_std = Sym_S2_Corr_HF_mean/Sym_S2_Corr_HF_std;
Sym_S1_Corr_SF_mean_i_std = Sym_S1_Corr_SF_mean/Sym_S1_Corr_SF_std;
Sym_S2_Corr_SF_mean_i_std = Sym_S2_Corr_SF_mean/Sym_S2_Corr_SF_std;
Sym_S1_Corr_UF_mean_i_std = Sym_S1_Corr_UF_mean/Sym_S1_Corr_UF_std;
Sym_S2_Corr_UF_mean_i_std = Sym_S2_Corr_UF_mean/Sym_S2_Corr_UF_std;

Sym_S1_Corr_MF_std_i_Sym_S2_Corr_MF_mean = Sym_S1_Corr_MF_mean_i_std/Sym_S2_Corr_MF_mean;
Sym_S1_Corr_mean_i_Sym_S2_Corr_mean=Sym_S1_Corr_mean/Sym_S2_Corr_mean;
Sym_S1_Corr_std_i_Sym_S2_Corr_std = Sym_S1_Corr_std/Sym_S2_Corr_std;

Sym_S1_std_mean = mean(symmetryS1std);
Sym_S2_std_mean = mean(symmetryS2std);

Sym_S1_std_std=std(symmetryS1std);
Sym_S2_std_std=std(symmetryS2std);

Skew_Glob_Raw = skewness(PCG_resampled);
Skew_Glob_LF = skewness(lowFEnvelope);
Skew_Glob_MF = skewness(midFEnvelope);
Skew_Glob_HF = skewness(highFEnvelope);
Skew_Glob_SF = skewness(superFEnvelope);
Skew_Glob_UF = skewness(ultraFEnvelope);

Skew_Loc_Raw_mean = mean(skewnessRAW);
Skew_Loc_Raw_std = std(skewnessRAW);

Skew_Loc_LF_mean = mean(skewnessLOW);
Skew_Loc_LF_std = std(skewnessLOW);

Skew_Loc_MF_std = std(skewnessMID);
Skew_Loc_MF_mean= mean(skewnessMID);

Skew_Loc_HF_mean = mean(skewnessHIGH);
Skew_Loc_HF_std = std(skewnessHIGH);

Skew_Loc_SF_mean = mean(skewnessSUPER);
Skew_Loc_SF_std = std(skewnessSUPER);

Skew_Loc_UF_mean = mean(skewnessULTRA);
Skew_Loc_UF_std = std(skewnessULTRA);

Skew_S1_LF_mean = mean(skewS1LOW);
Skew_S1_LF_std = std(skewS1LOW);
Skew_S2_LF_mean = mean(skewS2LOW);
Skew_S2_LF_std = std(skewS2LOW);

Skew_S1_MF_mean = mean(skewS1MID);
Skew_S1_MF_std = std(skewS1MID);
Skew_S2_MF_mean = mean(skewS2MID);
Skew_S2_MF_std = std(skewS2MID);

Skew_S1_HF_mean = mean(skewS1HIGH);
Skew_S1_HF_std = std(skewS1HIGH);
Skew_S2_HF_mean = mean(skewS2HIGH);
Skew_S2_HF_std = std(skewS2HIGH);

Skew_S1_SF_mean = mean(skewS1SUPER);
Skew_S1_SF_std = std(skewS1SUPER);
Skew_S2_SF_mean = mean(skewS2SUPER);
Skew_S2_SF_std = std(skewS2SUPER);

Skew_S1_UF_mean = mean(skewS1ULTRA);
Skew_S1_UF_std = std(skewS1ULTRA);
Skew_S2_UF_mean = mean(skewS2ULTRA);
Skew_S2_UF_std = std(skewS2ULTRA);

centerCorr = fix(length(avgLShape)/2);
offCorr = fix(centerCorr/2);
sc = centerCorr-offCorr;
ec = centerCorr+offCorr;

shortL = avgLShape(sc:ec);
shortM = avgMShape(sc:ec);
shortH = avgHShape(sc:ec);
shortS = avgSShape(sc:ec);
shortU = avgUShape(sc:ec);

shortX = (1:length(shortL))';

XT_L = mean(shortL.*shortX);
XT_M = mean(shortM.*shortX);
XT_H = mean(shortH.*shortX);
XT_S = mean(shortS.*shortX);
XT_U = mean(shortU.*shortX);