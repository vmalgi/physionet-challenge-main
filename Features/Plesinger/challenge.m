function classifyResult = challenge(recordName)
%
% Entry for the 2016 PhysioNet/CinC Challenge.
%
% F.Plesinger 2016

%% Load data and resample data
[PCG, Fs1, nbits1] = wavread([recordName '.wav']);  % load data


ultraFEnvelope = hobalka(PCG,Fs1,400,800);

Fs = 1000;
PCG_resampled      = resample(PCG,Fs,Fs1); % resample to schmidt_options.audio_Fs (1000 Hz)

ultraFEnvelope = resample(ultraFEnvelope,Fs,Fs1);

lfMin = 15; %15
lfMax = 90; %45

mfMin = 55;
mfMax = 150;

hfMin = 100;
hfMax = 250;

sfMin = 200;
sfMax = 450;

invalid_window=1000;
invalid_step = 500;
invalid = zeros(length(PCG_resampled),1);

% globalSTD = std(PCG_resampled);

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
env_60_200 = hobalka(PCG_resampled,Fs,60,200); %60-200
frequencyRatio = sum(env_60_200)/sum(env_1_40);

%stdev
% winSize = 50;
% % winStep = 10;
% stdev = zeros(length(PCG_resampled):1);
% for i=1:winSize:length(PCG_resampled)-winSize;
%    stdev(i:i+winSize)=std(PCG_resampled(i:i+winSize));
% end




% threshold = prctile(lowFEnvelope.75);

% samples = 1:length(PCG_resampled);

[s1positions s2positions extrema locExtrema minHistoVal maxHistoVal histogram avgLShape avgMShape avgHShape avgSShape avgUShape mdetcorrs fftcorrs] = getS1S2(PCG_resampled,Fs,lowFEnvelope,midFEnvelope,highFEnvelope,superFEnvelope, ultraFEnvelope, invalid );

rrs = diff(locExtrema(s1positions));

s1s2Widths = locExtrema(s2positions)-locExtrema(s1positions); %pro statistiku

expectedNumberOfBeats = length(PCG_resampled)/1000; %60 za minutu je pøedpokládaný poèet tepù. To je potøeba zpøesnit!
beatsCountRatio = length(s1positions)/expectedNumberOfBeats;

wrongForFeatureExtraction = (length(rrs)<5 | mean(rrs)<2*std(rrs) | mean(s1s2Widths)>600 | mean(rrs)>1600 | mean(rrs)<500 | median(rrs)<2*median(s1s2Widths));

%tøída souboru
class=0;

%use class - (training - a-e) differentiation:
if diffGroups

[pathstr,name,ext] = fileparts(recordName) ;

if strfind(name(1), 'a')
class=1;
end

if strfind(name(1), 'b')
class=2;
end

if strfind(name(1), 'c')
class=3;
end

if strfind(name(1), 'd')
class=4;
end

if strfind(name(1), 'e')
class=5;
end

if strfind(name(1), 'f')
class=6;
end

end

%oblasti v A.B(S1).C.D.E(S2).F.G - Extract features
ri=0;

esum = 0.6;
SK2M = 1.7;

skore = 0;

for index=1:length(s1positions)
   
   index = fix(index);
    
   s1pos = locExtrema(s1positions(index));
   s2pos = locExtrema(s2positions(index));
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
   
   if startF>length(PCG_resampled) | endF>length(PCG_resampled)
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

if doTable|doPictures
 %naètení odpovìdí
    [pathstr,flnm,ext] = fileparts(recordName);
    refFileName = [pathstr '\REFERENCE_new.CSV'];
    ftoread = refFileName;
    fid = fopen(ftoread);
    refValues = textscan(fid, '%s %f %f', 'Delimiter',','); 
    fclose (fid);
    numbers = refValues{2};
    sqi = refValues{3};
    thisIndex = find(strcmp(refValues{1},flnm));
    correctAnswer = refValues{2}(thisIndex);
    quality =  refValues{3}(thisIndex);
  % konec naètení správné odpovìdi
  
  %naètení patologie
   refFileName = [pathstr '\diag.csv'];
    ftoread = refFileName;
    fid = fopen(ftoread);
    refValues = textscan(fid, '%s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s %s', 'Delimiter',';'); 
    fclose (fid);
    thisIndex = find(strcmp(refValues{1},flnm));
    diagnoses = refValues{4};
    thisDiagnose =refValues{4}(thisIndex);
    
    diagnoseCode = getDiagnoseCode(thisDiagnose);
    
    
  %konec naètení patologie
  
end

skore=0;

Corr_FFT_mean = mean(mdetcorrs);
Corr_FFT_std = std(mdetcorrs);
Corr_FFT_skw = skewness(mdetcorrs);
Corr_FFT_kurt = kurtosis(mdetcorrs);
Corr_FFT_median = median(mdetcorrs);

Corr_FFT_mean_i_std = Corr_FFT_mean/Corr_FFT_std;



if exist('dataAmed')

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
%%%%%%

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
%%%%%%


if (minHistoVal~=maxHistoVal) & (length(s1positions)>2)

%výpoèet skore

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

else
    skore = 0;
end
end

classifyResult = 0;

limitCheck = 0;


if exist('Sym_S1_Corr_mean')
    
    %limit
    %check%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if SDMiEM>1.799603 || DMiFM>2.619169 || SDMiFM>3.842556 || SAMiCM<0.1964864 || DMiEM>1.254161
      classifyResult = 1;
      limitCheck = 1;
    end
    
    if classifyResult==0
    if Skew_Loc_UF_mean<0.6244103 || Skew_Glob_UF<0.6534703 || Corr_FFT_std<10.72538 || Skew_S2_LF_std<0.1060156 || UBMiBM<0.005140321
      classifyResult = -1;
      limitCheck = -1;
    end
    end
    
    
    if class==2 && classifyResult == 0
        if Sym_S2_Corr_SF_mean>0.5117335 || Sym_S2_Corr_MF_mean_i_std>5.642618 || DMiEM<0.1972007
           classifyResult = -1;
           limitCheck = -1; 
        end
    end
    
%     if class==4 && classifyResult ==0
%         if Corr_FFT_kurt>66.13723 || BMiCM<1.609617
%             classifyResult = -1;
%              limitCheck = -1;
%         end
%     end
    
    if class==5 && classifyResult == 0
        
          if UACDFMiBEM>2.428094 || UBMiCM<0.763442
              classifyResult = 1;
              limitCheck = 1;
         end
        
         if classifyResult == 0 && UDMiDM<0.073 || Skew_S1_UF_mean<0.6
             classifyResult = -1;
             limitCheck = -1;
         end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    if class==1 && classifyResult==0
%%%%%%% Autogenerated code from Histofind - 5.8.2016 - 14:42:59 %%%%%%%%
L1=0.19;L2=0.4299999;
GSKORE=0;
%SystoleWidth_mean
SystoleWidth_mean_PROB=[0.5,1,0,0,0.4,0.24,0.2727273,0.375,0.3783784,0.3508772,0.4545455,0.2698413,0.25,0.09677419,0.07692308,0,0,0,1,1];
SystoleWidth_mean_WEIGHT=[0.004889976,0.002444988,0.01466993,0.01711491,0.02444988,0.06112469,0.05378973,0.07823961,0.09046455,0.1393643,0.1075795,0.1540342,0.09779951,0.07579462,0.03178484,0.02689487,0.007334963,0.007334963,0.002444988,0.002444988];
DIVIDER=(428.9091-193.4342)/19;
INDEX=floor((SystoleWidth_mean-193.4342)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SystoleWidth_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SystoleWidth_mean_WEIGHT(INDEX);
end
end

%FM_mean
FM_mean_PROB=[0.2631579,0.7727273,NaN,0,0,NaN,NaN,NaN,0,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,0,0];
FM_mean_WEIGHT=[0.9290953,0.05378973,0,0.002444988,0.002444988,0,0,0,0.002444988,0,0,0,0.004889976,0,0,0,0,0,0.002444988,0.002444988];
DIVIDER=(0.0955-0.0003)/19;
INDEX=floor((FM_mean-0.0003)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =FM_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*FM_mean_WEIGHT(INDEX);
end
end

%AMiCM
AMiCM_PROB=[0,0,0,0,0,0.1333333,0.15,0.2352941,0.2432432,0.2903226,0.3877551,0.3793103,0.326087,0.3137255,0.2380952,0.375,0.3636364,0.4,0.3333333,NaN];
AMiCM_WEIGHT=[0.007334963,0.01711491,0.01466993,0.01711491,0.01466993,0.03667482,0.04889975,0.04156479,0.09046455,0.07579462,0.1198044,0.1418093,0.1124694,0.1246944,0.05134474,0.03911981,0.02689487,0.01222494,0.007334963,0];
DIVIDER=(1.3709-0.0488)/19;
INDEX=floor((AMiCM-0.0488)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =AMiCM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*AMiCM_WEIGHT(INDEX);
end
end

%DMiFM
DMiFM_PROB=[0.3301887,0.3519313,0,0,0,0,0,0,0,NaN,0,0,0,0,0,0,NaN,0,0,NaN];
DMiFM_WEIGHT=[0.2591687,0.5696822,0.06112469,0.03667482,0.007334963,0.01711491,0.009779952,0.01466993,0.002444988,0,0.002444988,0.002444988,0.002444988,0.002444988,0.002444988,0.002444988,0,0.004889976,0.002444988,0];
DIVIDER=(7.7547-0.5805)/19;
INDEX=floor((DMiFM-0.5805)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =DMiFM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*DMiFM_WEIGHT(INDEX);
end
end

%AMiGM
AMiGM_PROB=[0,NaN,0,0,0,0.2,0.2875,0.2992701,0.3239437,0.2777778,0.2413793,0.8333333,0,0.6666667,1,0,NaN,NaN,0,0];
AMiGM_WEIGHT=[0.002444988,0,0.002444988,0.009779952,0.02200489,0.06112469,0.195599,0.3349633,0.1735941,0.08801956,0.07090464,0.01466993,0.007334963,0.007334963,0.002444988,0.002444988,0,0,0.002444988,0.002444988];
DIVIDER=(2.4116-0.2601)/19;
INDEX=floor((AMiGM-0.2601)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =AMiGM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*AMiGM_WEIGHT(INDEX);
end
end

%UAMiBM
UAMiBM_PROB=[0,0.4,0.4166667,0.3888889,0.4057971,0.2454545,0.2464789,0.1304348,0.1666667,0,1,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1];
UAMiBM_WEIGHT=[0.007334963,0.01222494,0.02933985,0.08801956,0.1687042,0.2689486,0.3471883,0.05623472,0.01466993,0.002444988,0.002444988,0,0,0,0,0,0,0,0,0.002444988];
DIVIDER=(2.7026-0.0678)/19;
INDEX=floor((UAMiBM-0.0678)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =UAMiBM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*UAMiBM_WEIGHT(INDEX);
end
end

%Sym_S2_Corr_MF_std
Sym_S2_Corr_MF_std_PROB=[0.2,0,0.2666667,0.4736842,0.4583333,0.3043478,0.3448276,0.3888889,0.254902,0.2444444,0.2857143,0.2926829,0.2083333,0.0952381,0.1538462,0.1666667,1,NaN,0,NaN];
Sym_S2_Corr_MF_std_WEIGHT=[0.01222494,0.01466993,0.03667482,0.04645477,0.05867971,0.05623472,0.07090464,0.08801956,0.1246944,0.1100245,0.1198044,0.1002445,0.05867971,0.05134474,0.03178484,0.01466993,0.002444988,0,0.002444988,0];
DIVIDER=(0.5269-0.0525)/19;
INDEX=floor((Sym_S2_Corr_MF_std-0.0525)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Sym_S2_Corr_MF_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Sym_S2_Corr_MF_std_WEIGHT(INDEX);
end
end

%Sym_S2_Corr_SF_std
Sym_S2_Corr_SF_std_PROB=[0,0.2,0.3636364,0.25,0.25,0.2083333,0.2173913,0.25,0.3409091,0.26,0.4186046,0.4102564,0.25,0.3478261,0.08333334,0.2857143,0,0,0,NaN];
Sym_S2_Corr_SF_std_WEIGHT=[0.004889976,0.01222494,0.02689487,0.02933985,0.07823961,0.05867971,0.05623472,0.07823961,0.1075795,0.1222494,0.1051345,0.09535452,0.07823961,0.05623472,0.02933985,0.03422983,0.01466993,0.002444988,0.009779952,0];
DIVIDER=(0.3777-0.1231)/19;
INDEX=floor((Sym_S2_Corr_SF_std-0.1231)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Sym_S2_Corr_SF_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Sym_S2_Corr_SF_std_WEIGHT(INDEX);
end
end

%Skew_Glob_MF
Skew_Glob_MF_PROB=[0.1666667,0.2535211,0.3033708,0.31,0.44,0,0.1,0,0,0,0,1,NaN,0,NaN,NaN,NaN,0,0,NaN];
Skew_Glob_MF_WEIGHT=[0.01466993,0.1735941,0.4352078,0.2444988,0.06112469,0.01466993,0.02444988,0.007334963,0.002444988,0.002444988,0.009779952,0.002444988,0,0.002444988,0,0,0,0.002444988,0.002444988,0];
DIVIDER=(22.6563-0.9273)/19;
INDEX=floor((Skew_Glob_MF-0.9273)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Glob_MF_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Glob_MF_WEIGHT(INDEX);
end
end

%Skew_Loc_SF_std
Skew_Loc_SF_std_PROB=[0.1666667,0.4285714,0.3,0.2258064,0.2083333,0.2962963,0.372549,0.2181818,0.1956522,0.28125,0.4,0.2307692,0.5454546,0.4,0.6,0,0,NaN,1,1];
Skew_Loc_SF_std_WEIGHT=[0.01466993,0.01711491,0.02444988,0.07579462,0.1173594,0.1320293,0.1246944,0.1344743,0.1124694,0.07823961,0.07334963,0.03178484,0.02689487,0.01222494,0.01222494,0.004889976,0.002444988,0,0.002444988,0.002444988];
DIVIDER=(3.2772-0.2352)/19;
INDEX=floor((Skew_Loc_SF_std-0.2352)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Loc_SF_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Loc_SF_std_WEIGHT(INDEX);
end
end

%Skew_Loc_UF_mean
Skew_Loc_UF_mean_PROB=[0,0,0.1,0,0.07407407,0.07692308,0.28125,0.2727273,0.295082,0.4363636,0.4210526,0.3947369,0.3333333,0.6,0.25,1,1,0,1,NaN];
Skew_Loc_UF_mean_WEIGHT=[0.007334963,0.01711491,0.02444988,0.06112469,0.06601467,0.06356968,0.07823961,0.1075795,0.1491442,0.1344743,0.09290954,0.09290954,0.05134474,0.01222494,0.0195599,0.009779952,0.002444988,0.007334963,0.002444988,0];
DIVIDER=(9.3054-0.5453)/19;
INDEX=floor((Skew_Loc_UF_mean-0.5453)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Loc_UF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Loc_UF_mean_WEIGHT(INDEX);
end
end

%Skew_S1_LF_std
Skew_S1_LF_std_PROB=[0.4,0.5384616,0.2272727,0.5142857,0.3714286,0.3559322,0.2222222,0.2291667,0.1764706,0.3076923,0.3333333,0.09090909,0.08333334,0.1,0.2,NaN,NaN,1,0,NaN];
Skew_S1_LF_std_WEIGHT=[0.01222494,0.03178484,0.05378973,0.08557457,0.08557457,0.1442543,0.1540342,0.1173594,0.08312958,0.06356968,0.04400978,0.05378973,0.02933985,0.02444988,0.01222494,0,0,0.002444988,0.002444988,0];
DIVIDER=(0.8975-0.078)/19;
INDEX=floor((Skew_S1_LF_std-0.078)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S1_LF_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S1_LF_std_WEIGHT(INDEX);
end
end

%Skew_S2_LF_mean
Skew_S2_LF_mean_PROB=[0,NaN,0,0.1818182,0,0.4,0.2,0.2083333,0.3043478,0.2439024,0.3043478,0.4313726,0.3953488,0.25,0.4411765,0.1538462,0.1875,0.1666667,0.25,0];
Skew_S2_LF_mean_WEIGHT=[0.002444988,0,0.01222494,0.02689487,0.03911981,0.02444988,0.03667482,0.05867971,0.05623472,0.1002445,0.1124694,0.1246944,0.1051345,0.08801956,0.08312958,0.06356968,0.03911981,0.01466993,0.009779952,0.002444988];
DIVIDER=(2.9139-0.3244)/19;
INDEX=floor((Skew_S2_LF_mean-0.3244)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S2_LF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S2_LF_mean_WEIGHT(INDEX);
end
end

%Skew_S2_LF_std
Skew_S2_LF_std_PROB=[0.25,0.4285714,0.4375,0.5263158,0.2903226,0.4444444,0.475,0.12,0.2380952,0.3170732,0.1,0.2580645,0.1363636,0.3157895,0.5,0,0,0,0,NaN];
Skew_S2_LF_std_WEIGHT=[0.009779952,0.01711491,0.03911981,0.04645477,0.07579462,0.08801956,0.09779951,0.1222494,0.1026895,0.1002445,0.09779951,0.07579462,0.05378973,0.04645477,0.009779952,0.004889976,0.004889976,0.002444988,0.004889976,0];
DIVIDER=(0.8743-0.0854)/19;
INDEX=floor((Skew_S2_LF_std-0.0854)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S2_LF_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S2_LF_std_WEIGHT(INDEX);
end
end

%XT_M
XT_M_PROB=[0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,1,1,0.2831633,0.1666667,NaN,1,1];
XT_M_WEIGHT=[0.002444988,0,0,0,0,0,0,0,0,0,0,0,0,0.002444988,0.002444988,0.9584352,0.02933985,0,0.002444988,0.002444988];
DIVIDER=(987315.4--4823000)/19;
INDEX=floor((XT_M--4823000)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =XT_M_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*XT_M_WEIGHT(INDEX);
end
end

if GSKORE>-0.0005593141
classifyResult=-1;
end

if GSKORE<-0.0005593141
classifyResult=1;
end

%%%%%%%%% END AUTOGENERATED CODE FROM Histofind %%%%%%%%%%%%%%%%
    end

    
    if class==2 && classifyResult==0
%%%%%%% Autogenerated code from Histofind - 5.8.2016 - 10:30:10 %%%%%%%%
L1=0.6799996;L2=0.8499995;
GSKORE=0;
%SystoleWidth_mean
SystoleWidth_mean_PROB=[1,0.4,0.7272727,0.7142857,0.72,0.72,0.88,0.8064516,0.7121212,0.8055556,0.9090909,0.8421053,0.7586207,0.8823529,0.5,1,0,NaN,NaN,1];
SystoleWidth_mean_WEIGHT=[0.004081632,0.01020408,0.02244898,0.02857143,0.05102041,0.05102041,0.1020408,0.1265306,0.1346939,0.1469388,0.1122449,0.07755102,0.05918367,0.03469388,0.02857143,0.004081632,0.004081632,0,0,0.002040816];
DIVIDER=(445.4-202.7273)/19;
INDEX=floor((SystoleWidth_mean-202.7273)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SystoleWidth_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SystoleWidth_mean_WEIGHT(INDEX);
end
end

%SystoleWidth_std
SystoleWidth_std_PROB=[0.7765958,0.875,0.9245283,0.8,0.7619048,0.7631579,0.7567568,0.7333333,0.7222222,0.6666667,0.7222222,0.6666667,0.4,0.6666667,0,1,NaN,1,NaN,1];
SystoleWidth_std_WEIGHT=[0.1918367,0.1469388,0.1081633,0.09183674,0.08571429,0.07755102,0.0755102,0.06122449,0.03673469,0.04897959,0.03673469,0.0122449,0.01020408,0.006122449,0.002040816,0.002040816,0,0.004081632,0,0.002040816];
DIVIDER=(167.5179-2.1213)/19;
INDEX=floor((SystoleWidth_std-2.1213)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SystoleWidth_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SystoleWidth_std_WEIGHT(INDEX);
end
end

%GS_mean
GS_mean_PROB=[0.8571429,0.7954546,0.8688524,0.75,0.744186,0.7027027,0.5714286,0.75,1,0.5,1,1,1,0.6666667,0.5,NaN,0.5,NaN,NaN,1];
GS_mean_WEIGHT=[0.05714286,0.1795918,0.2489796,0.1795918,0.0877551,0.0755102,0.05714286,0.03265306,0.02040816,0.0122449,0.01428571,0.008163265,0.01020408,0.006122449,0.004081632,0,0.004081632,0,0,0.002040816];
DIVIDER=(3.3321-0.0518)/19;
INDEX=floor((GS_mean-0.0518)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =GS_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*GS_mean_WEIGHT(INDEX);
end
end

%BD_mean
BD_mean_PROB=[0.8769231,0.8210526,0.8041237,0.6206896,0.7045454,0.8666667,0.7391304,0.8888889,0.8571429,0.8571429,0.75,0.6666667,1,0.4,1,NaN,NaN,NaN,0,NaN];
BD_mean_WEIGHT=[0.1326531,0.1938775,0.1979592,0.1183673,0.08979592,0.09183674,0.04693878,0.03673469,0.02857143,0.01428571,0.01632653,0.0122449,0.004081632,0.01020408,0.004081632,0,0,0,0.002040816,0];
DIVIDER=(0.0562-0.0009)/19;
INDEX=floor((BD_mean-0.0009)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =BD_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*BD_mean_WEIGHT(INDEX);
end
end

%ADiGD
ADiGD_PROB=[0.7037037,0.7525773,0.7709923,0.7234042,0.8867925,0.7878788,0.9047619,0.9230769,1,1,1,1,1,1,1,NaN,NaN,NaN,NaN,1];
ADiGD_WEIGHT=[0.05510204,0.1979592,0.2673469,0.1918367,0.1081633,0.06734694,0.04285714,0.02653061,0.01428571,0.006122449,0.008163265,0.002040816,0.002040816,0.002040816,0.006122449,0,0,0,0,0.002040816];
DIVIDER=(6.0132-0.113)/19;
INDEX=floor((ADiGD-0.113)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =ADiGD_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*ADiGD_WEIGHT(INDEX);
end
end

%AMiCM
AMiCM_PROB=[0,0.5,0.4615385,0.4166667,0.7857143,0.8,0.7818182,0.8695652,0.8888889,0.9142857,0.8571429,0.8888889,1,1,1,NaN,1,NaN,0,1];
AMiCM_WEIGHT=[0.01632653,0.008163265,0.02653061,0.04897959,0.08571429,0.1632653,0.2244898,0.1408163,0.1469388,0.07142857,0.02857143,0.01836735,0.004081632,0.004081632,0.004081632,0,0.004081632,0,0.002040816,0.002040816];
DIVIDER=(2.0253-0.2164)/19;
INDEX=floor((AMiCM-0.2164)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =AMiCM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*AMiCM_WEIGHT(INDEX);
end
end

%UAS_mean
UAS_mean_PROB=[0.7735849,0.8083333,0.8181818,0.7586207,0.8421053,0.7647059,1,1,0.6666667,0.5,1,0,1,1,NaN,0,NaN,NaN,1,NaN];
UAS_mean_WEIGHT=[0.4326531,0.244898,0.1346939,0.05918367,0.03877551,0.03469388,0.01020408,0.01020408,0.006122449,0.004081632,0.004081632,0.008163265,0.006122449,0.002040816,0,0.002040816,0,0,0.002040816,0];
DIVIDER=(0.6458-0.0045)/19;
INDEX=floor((UAS_mean-0.0045)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =UAS_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*UAS_mean_WEIGHT(INDEX);
end
end

%UCS_mean
UCS_mean_PROB=[0.8175182,0.8296296,0.7692308,0.7619048,0.7714286,0.8,0.4375,0.75,0.7777778,1,1,0.5,0.6666667,NaN,0.5,NaN,1,NaN,0.5,NaN];
UCS_mean_WEIGHT=[0.2795918,0.2755102,0.1591837,0.08571429,0.07142857,0.02040816,0.03265306,0.0244898,0.01836735,0.004081632,0.006122449,0.004081632,0.006122449,0,0.004081632,0,0.004081632,0,0.004081632,0];
DIVIDER=(0.5355-0.0036)/19;
INDEX=floor((UCS_mean-0.0036)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =UCS_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*UCS_mean_WEIGHT(INDEX);
end
end

%Sym_S1_Corr_mean
Sym_S1_Corr_mean_PROB=[1,0.5,1,1,0.6,0.4,0.6,0.8235294,0.84,0.78125,0.8648649,0.8636364,0.7636364,0.8,0.8,0.8644068,0.6341463,0.7804878,0.7222222,1];
Sym_S1_Corr_mean_WEIGHT=[0.002040816,0.004081632,0.002040816,0.002040816,0.01020408,0.01020408,0.02040816,0.03469388,0.05102041,0.06530612,0.0755102,0.08979592,0.1122449,0.1020408,0.09183674,0.1204082,0.08367347,0.08367347,0.03673469,0.002040816];
DIVIDER=(0.9778-0.2577)/19;
INDEX=floor((Sym_S1_Corr_mean-0.2577)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Sym_S1_Corr_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Sym_S1_Corr_mean_WEIGHT(INDEX);
end
end

%Sym_S2_Corr_mean
Sym_S2_Corr_mean_PROB=[1,NaN,1,NaN,NaN,NaN,0.75,0.3333333,0.25,0.5454546,0.8333333,0.7692308,0.7291667,0.7,0.7605634,0.8550724,0.8533334,0.90625,0.8333333,NaN];
Sym_S2_Corr_mean_WEIGHT=[0.002040816,0,0.002040816,0,0,0,0.008163265,0.006122449,0.01632653,0.02244898,0.0244898,0.07959183,0.09795918,0.122449,0.144898,0.1408163,0.1530612,0.1306122,0.04897959,0];
DIVIDER=(0.9735--0.0157)/19;
INDEX=floor((Sym_S2_Corr_mean--0.0157)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Sym_S2_Corr_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Sym_S2_Corr_mean_WEIGHT(INDEX);
end
end

%Sym_S2_Corr_std
Sym_S2_Corr_std_PROB=[0.7333333,0.8823529,0.9117647,0.8571429,0.8717949,0.8095238,0.7884616,0.8292683,0.6888889,0.8148148,0.6428571,0.7083333,0.7916667,0.75,0.6923077,0.5,0.25,0.6,1,1];
Sym_S2_Corr_std_WEIGHT=[0.03061225,0.06938776,0.06938776,0.08571429,0.07959183,0.08571429,0.1061224,0.08367347,0.09183674,0.05510204,0.05714286,0.04897959,0.04897959,0.01632653,0.02653061,0.01632653,0.008163265,0.01020408,0.006122449,0.004081632];
DIVIDER=(0.4918-0.0114)/19;
INDEX=floor((Sym_S2_Corr_std-0.0114)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Sym_S2_Corr_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Sym_S2_Corr_std_WEIGHT(INDEX);
end
end

%Skew_Glob_SF
Skew_Glob_SF_PROB=[0.05882353,0.8611111,0.7808219,0.8596491,0.8367347,0.7619048,0.9,0.825,0.9032258,0.7647059,0.7,0.75,0.8333333,0.7647059,0.6666667,0.8333333,0.75,0.6666667,1,NaN];
Skew_Glob_SF_WEIGHT=[0.03469388,0.07346939,0.1489796,0.1163265,0.1,0.08571429,0.06122449,0.08163265,0.06326531,0.03469388,0.04081633,0.03265306,0.0244898,0.03469388,0.01836735,0.0122449,0.008163265,0.01836735,0.01020408,0];
DIVIDER=(31.8994-1.8875)/19;
INDEX=floor((Skew_Glob_SF-1.8875)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Glob_SF_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Glob_SF_WEIGHT(INDEX);
end
end

%Skew_Loc_MF_mean
Skew_Loc_MF_mean_PROB=[1,NaN,0.4545455,0.7,0.7272727,0.65625,0.5964912,0.9387755,0.8333333,0.8923077,0.8644068,0.8048781,0.7916667,0.9047619,0.7857143,0.8333333,0.1666667,0.6666667,1,0];
Skew_Loc_MF_mean_WEIGHT=[0.004081632,0,0.02244898,0.02040816,0.04489796,0.06530612,0.1163265,0.1,0.1346939,0.1326531,0.1204082,0.08367347,0.04897959,0.04285714,0.02857143,0.0122449,0.0122449,0.006122449,0.002040816,0.002040816];
DIVIDER=(4.3556-0.9733)/19;
INDEX=floor((Skew_Loc_MF_mean-0.9733)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Loc_MF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Loc_MF_mean_WEIGHT(INDEX);
end
end

%Skew_S2_LF_std
Skew_S2_LF_std_PROB=[0.7777778,0.8260869,0.8461539,0.8372093,0.82,0.7916667,0.7846154,0.8333333,0.84,0.8205128,0.5909091,0.8148148,0.4666667,0.8571429,0.8,0,0,NaN,NaN,1];
Skew_S2_LF_std_WEIGHT=[0.01836735,0.04693878,0.05306122,0.0877551,0.1020408,0.09795918,0.1326531,0.1102041,0.1020408,0.07959183,0.04489796,0.05510204,0.03061225,0.01428571,0.01020408,0.008163265,0.004081632,0,0,0.002040816];
DIVIDER=(1.1527-0.0938)/19;
INDEX=floor((Skew_S2_LF_std-0.0938)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S2_LF_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S2_LF_std_WEIGHT(INDEX);
end
end

%Skew_S2_MF_mean
Skew_S2_MF_mean_PROB=[0.7777778,0.4285714,0.4,0.5714286,0.6571429,0.6046512,0.7446808,0.8461539,0.8928571,0.8035714,0.8947368,0.9705882,0.8260869,0.8636364,0.8666667,1,1,1,1,NaN];
Skew_S2_MF_mean_WEIGHT=[0.01836735,0.01428571,0.03061225,0.04285714,0.07142857,0.0877551,0.09591836,0.1061224,0.1142857,0.1142857,0.07755102,0.06938776,0.04693878,0.04489796,0.03061225,0.008163265,0.0122449,0.0122449,0.002040816,0];
DIVIDER=(3.4113-0.8624)/19;
INDEX=floor((Skew_S2_MF_mean-0.8624)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S2_MF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S2_MF_mean_WEIGHT(INDEX);
end
end

%Skew_S2_SF_mean
Skew_S2_SF_mean_PROB=[1,0.25,0.2857143,0.6666667,0.6296296,0.6,0.7,0.796875,0.7941176,0.8809524,0.8333333,0.8666667,0.9655172,0.9310345,1,1,0.5,0.25,0.3333333,0];
Skew_S2_SF_mean_WEIGHT=[0.004081632,0.008163265,0.01428571,0.01836735,0.05510204,0.06122449,0.08163265,0.1306122,0.1387755,0.08571429,0.1346939,0.09183674,0.05918367,0.05918367,0.02040816,0.01632653,0.004081632,0.008163265,0.006122449,0.002040816];
DIVIDER=(5.4084-0.6858)/19;
INDEX=floor((Skew_S2_SF_mean-0.6858)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S2_SF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S2_SF_mean_WEIGHT(INDEX);
end
end

if GSKORE>0.0006947567
classifyResult=-1;
end

if GSKORE<0.0006947567
classifyResult=1;
end

%%%%%%%%% END AUTOGENERATED CODE FROM Histofind %%%%%%%%%%%%%%%%
    end
    
    if class==3  && classifyResult==0
        
        %%% ruènì urèená náhrada
        
        if Skew_S2_HF_mean > 2.46
            classifyResult = -1;
        else
            classifyResult = 1;
        end
        
        %%%%%%% Autogenerated code from STATFIND 26.7.2016 - 8:59:06 %%%%%%%%
% L1=0.11;L2=0.6499997;
% GSKORE=0;
% %Skew_S2_HF_mean
% Skew_S2_HF_mean_PROB=[0,0,0,0,0,0,0,0,0,0,NaN,1,1,1,1,1,1,1,1,1]; %ruènì upraveno!!!!!!!!
% Skew_S2_HF_mean_WEIGHT=[0.06666667,0.06666667,0.06666667,0.06666667,0.1,0.1,0,0.1,0.1,0.1,0.1333333,0,0.06666667,0.06666667,0.06666667,0.03333334,0.03333334,0.03333334,0.03333334,0.03333334,0.03333334]; %ruènì upraqveno!!!!
% DIVIDER=(3.8257-0.8405)/19;
% INDEX=floor((Skew_S2_HF_mean-0.8405)/DIVIDER)+1;
% if INDEX>0 & INDEX<20
% MP=0;
% PROB =Skew_S2_HF_mean_PROB(INDEX);
% if PROB<L1
% MP=-(L1-PROB);
% end
% if PROB>L2
% MP=PROB-L2;
% end
% if ~isnan(MP)
% GSKORE=GSKORE+MP*Skew_S2_HF_mean_WEIGHT(INDEX);
% end
% end
% 
% if GSKORE>0.0002353038
% classifyResult=-1;
% end
% 
% if GSKORE<0.0002353038
% classifyResult=1;
% end

%%%%%%%%% END AUTOGENERATED CODE FROM STATFIND %%%%%%%%%%%%%%%%
    end
    
    if class==4  && classifyResult==0
       %%%%%%% Autogenerated code from Histofind - 2.8.2016 - 10:47:01 %%%%%%%%
L1=0.14;L2=0.9999993;
GSKORE=0;
%BMiCM
BMiCM_PROB=[0,0,0,0.3333333,0.7777778,0.4,0.5,0.75,0,1,1,0.5,1,0,NaN,NaN,NaN,NaN,1,NaN];
BMiCM_WEIGHT=[0.07272727,0.09090909,0.03636364,0.05454545,0.1636364,0.09090909,0.1090909,0.07272727,0.07272727,0.07272727,0.05454545,0.03636364,0.01818182,0.01818182,0,0,0,0,0.03636364,0];
DIVIDER=(7.7288-0.2802)/19;
INDEX=floor((BMiCM-0.2802)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =BMiCM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*BMiCM_WEIGHT(INDEX);
end
end

%SAMiBM
SAMiBM_PROB=[0.5,1,0.625,1,0,0.3333333,0,0,1,0,0,0.5,0.5,0,1,0.6666667,NaN,0,1,1];
SAMiBM_WEIGHT=[0.03636364,0.01818182,0.1454545,0.07272727,0.03636364,0.05454545,0.07272727,0.05454545,0.03636364,0.05454545,0.05454545,0.07272727,0.07272727,0.01818182,0.09090909,0.05454545,0,0.01818182,0.01818182,0.01818182];
DIVIDER=(1.1257-0.2293)/19;
INDEX=floor((SAMiBM-0.2293)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SAMiBM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SAMiBM_WEIGHT(INDEX);
end
end

%Sym_S2_Sum_std
Sym_S2_Sum_std_PROB=[0,0.7,0.8,0.375,0.6,1,0,0,0.5,0.5,0,1,NaN,1,NaN,NaN,NaN,NaN,NaN,1];
Sym_S2_Sum_std_WEIGHT=[0.1272727,0.1818182,0.09090909,0.1454545,0.09090909,0.07272727,0.05454545,0.07272727,0.03636364,0.03636364,0.01818182,0.03636364,0,0.01818182,0,0,0,0,0,0.01818182];
DIVIDER=(1.0902-0.0589)/19;
INDEX=floor((Sym_S2_Sum_std-0.0589)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Sym_S2_Sum_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Sym_S2_Sum_std_WEIGHT(INDEX);
end
end

if GSKORE>-0.001490203
classifyResult=-1;
end

if GSKORE<-0.001490203
classifyResult=1;
end

%%%%%%%%% END AUTOGENERATED CODE FROM Histofind %%%%%%%%%%%%%%%%
    end
    
    if class==5 && classifyResult==0
        
%%%%%%% Autogenerated code from Histofind - 5.8.2016 - 15:42:08 %%%%%%%%
L1=0.6499997;L2=0.9999993;
GSKORE=0;
%SystoleWidth_std
SystoleWidth_std_PROB=[0.9929078,0.976923,0.8980892,0.8421053,0.6,0.48,0.25,0.4,0.1428571,0.6,0,0.5,NaN,NaN,NaN,1,NaN,NaN,NaN,1];
SystoleWidth_std_WEIGHT=[0.5038285,0.331802,0.08014293,0.02909648,0.02041858,0.01276161,0.008167433,0.005104645,0.003573252,0.002552323,0.0005104645,0.001020929,0,0,0,0.0005104645,0,0,0,0.0005104645];
DIVIDER=(134.9035-0.6467)/19;
INDEX=floor((SystoleWidth_std-0.6467)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SystoleWidth_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SystoleWidth_std_WEIGHT(INDEX);
end
end

%SDDiED
SDDiED_PROB=[0.9945652,0.9734849,0.9812734,0.9766537,0.970297,0.9698276,0.9509804,0.9236111,0.78,0.6842105,0.1333333,0.125,0.1428571,NaN,0,0,NaN,0,NaN,0];
SDDiED_WEIGHT=[0.09392547,0.1347626,0.136294,0.1311894,0.1546707,0.1184278,0.1041348,0.07350689,0.02552323,0.009698826,0.007656968,0.004083716,0.003573252,0,0.001020929,0.0005104645,0,0.0005104645,0,0.0005104645];
DIVIDER=(2.2709-0.0552)/19;
INDEX=floor((SDDiED-0.0552)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SDDiED_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SDDiED_WEIGHT(INDEX);
end
end

%SAMiBM
SAMiBM_PROB=[1,0.9230769,0.9032258,0.902439,0.984375,0.9649123,0.9626556,0.9900498,0.9778481,0.9370629,0.8333333,0.1176471,0.1428571,0,0,0,0,0,0,NaN];
SAMiBM_WEIGHT=[0.003062787,0.006636039,0.0158244,0.02092905,0.03266973,0.05819296,0.123022,0.3078101,0.3226136,0.07299643,0.009188361,0.008677897,0.003573252,0.004594181,0.003573252,0.003573252,0.0005104645,0.002041858,0.0005104645,0];
DIVIDER=(2.1801-0.0444)/19;
INDEX=floor((SAMiBM-0.0444)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SAMiBM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SAMiBM_WEIGHT(INDEX);
end
end

%UAMiGM
UAMiGM_PROB=[1,0.5,0.3333333,0.6,0.7234042,0.9337979,0.9790146,0.957346,0.8113208,0.5384616,0.75,0,1,1,0,1,NaN,0,1,0];
UAMiGM_WEIGHT=[0.001020929,0.001020929,0.001531394,0.007656968,0.02399183,0.1465033,0.5594691,0.215416,0.02705462,0.006636039,0.004083716,0.0005104645,0.001531394,0.001020929,0.0005104645,0.0005104645,0,0.0005104645,0.0005104645,0.0005104645];
DIVIDER=(1.7414-0.6157)/19;
INDEX=floor((UAMiGM-0.6157)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =UAMiGM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*UAMiGM_WEIGHT(INDEX);
end
end

%UACDFMiBEM
UACDFMiBEM_PROB=[0.5,1,1,1,0.8571429,0,0.7692308,1,0.9589041,0.9979101,0.9899713,0.6875,0.2307692,0,0,0,0,0,0,NaN];
UACDFMiBEM_WEIGHT=[0.001020929,0.001531394,0.001020929,0.001531394,0.003573252,0.0005104645,0.006636039,0.009698826,0.07452782,0.4885145,0.3563042,0.01633487,0.006636039,0.008677897,0.006636039,0.005104645,0.004594181,0.004083716,0.003062787,0];
DIVIDER=(3.2196-0.5913)/19;
INDEX=floor((UACDFMiBEM-0.5913)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =UACDFMiBEM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*UACDFMiBEM_WEIGHT(INDEX);
end
end

%Skew_Glob_Raw
Skew_Glob_Raw_PROB=[0,NaN,NaN,0,0,0.5,NaN,0,NaN,0.3333333,NaN,0.8,1,0.9565217,0.9596542,0.9447078,0.9610389,0.9705882,1,NaN];
Skew_Glob_Raw_WEIGHT=[0.0005104645,0,0,0.0005104645,0.0005104645,0.001020929,0,0.0005104645,0,0.001531394,0,0.002552323,0.006125574,0.02348137,0.1771312,0.6462481,0.1179173,0.01735579,0.004594181,0];
DIVIDER=(4.3017--19.7438)/19;
INDEX=floor((Skew_Glob_Raw--19.7438)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Glob_Raw_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Glob_Raw_WEIGHT(INDEX);
end
end

%Skew_Loc_UF_mean
Skew_Loc_UF_mean_PROB=[0.9969088,0.9373494,0.7777778,0.6585366,0.6,0.5384616,0.6666667,0.3333333,0.5555556,0.6666667,0.6666667,0,0,0.6,1,1,0.8571429,1,1,NaN];
Skew_Loc_UF_mean_WEIGHT=[0.6605411,0.2118428,0.05513017,0.02092905,0.01020929,0.006636039,0.004594181,0.003062787,0.004594181,0.003062787,0.003062787,0.001020929,0.001020929,0.002552323,0.003062787,0.001531394,0.003573252,0.002041858,0.001531394,0];
DIVIDER=(6.8111-0.3463)/19;
INDEX=floor((Skew_Loc_UF_mean-0.3463)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Loc_UF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Loc_UF_mean_WEIGHT(INDEX);
end
end

%Skew_S1_UF_mean
Skew_S1_UF_mean_PROB=[1,0.9915515,0.7935484,0.6875,0.7272727,0.55,0.5714286,0.5,0.8,0.25,0.625,0.5,0.6,NaN,0.5,1,0.6666667,1,0.6666667,1];
Skew_S1_UF_mean_WEIGHT=[0.16488,0.6646248,0.079122,0.03266973,0.01684533,0.01020929,0.003573252,0.005104645,0.002552323,0.002041858,0.004083716,0.002041858,0.002552323,0,0.002041858,0.002041858,0.001531394,0.002041858,0.001531394,0.0005104645];
DIVIDER=(5.003-0.2799)/19;
INDEX=floor((Skew_S1_UF_mean-0.2799)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S1_UF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S1_UF_mean_WEIGHT(INDEX);
end
end

%Skew_S1_UF_std
Skew_S1_UF_std_PROB=[0.9862638,0.9828031,0.9578313,0.8421053,0.8478261,0.7407407,0.8888889,0.6333333,0.3333333,0.8095238,0.5833333,0.5833333,0.6,0.6,0,0.3333333,NaN,NaN,NaN,1];
Skew_S1_UF_std_WEIGHT=[0.1858091,0.5936702,0.08473711,0.02909648,0.02348137,0.01378254,0.01378254,0.01531394,0.007656968,0.01071975,0.006125574,0.006125574,0.002552323,0.002552323,0.002552323,0.001531394,0,0,0,0.0005104645];
DIVIDER=(3.1161-0.0341)/19;
INDEX=floor((Skew_S1_UF_std-0.0341)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S1_UF_std_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S1_UF_std_WEIGHT(INDEX);
end
end

%Skew_S2_UF_mean
Skew_S2_UF_mean_PROB=[1,0.9891618,0.9074074,0.7578948,0.6744186,0.7333333,0.6666667,0.3333333,0.5,0.7142857,NaN,1,0,0.6666667,NaN,NaN,1,1,1,1];
Skew_S2_UF_mean_WEIGHT=[0.02756508,0.7064829,0.1653905,0.04849413,0.02194997,0.007656968,0.004594181,0.006125574,0.002041858,0.003573252,0,0.0005104645,0.001020929,0.001531394,0,0,0.0005104645,0.001531394,0.0005104645,0.0005104645];
DIVIDER=(4.6701-0.199)/19;
INDEX=floor((Skew_S2_UF_mean-0.199)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S2_UF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S2_UF_mean_WEIGHT(INDEX);
end
end

if GSKORE>-0.001395647
classifyResult=-1;
end

if GSKORE<-0.001395647
classifyResult=1;
end

%%%%%%%%% END AUTOGENERATED CODE FROM Histofind %%%%%%%%%%%%%%%%
    end
    
    if class==0  && classifyResult==0
  %%%%%%% Autogenerated code from Histofind - 5.8.2016 - 13:10:15 %%%%%%%%
L1=0.6699997;L2=0.9399994;
GSKORE=0;
%Sys_i_dia_mean
Sys_i_dia_mean_PROB=[0.827931,0.4144144,0.5833333,0.5,0.5,0.6666667,0,NaN,NaN,NaN,NaN,0,NaN,NaN,NaN,NaN,0,NaN,1,NaN];
Sys_i_dia_mean_WEIGHT=[0.9514436,0.03641732,0.003937008,0.00328084,0.002624672,0.000984252,0.000328084,0,0,0,0,0.000328084,0,0,0,0,0.000328084,0,0.000328084,0];
DIVIDER=(53.8283-0.3379)/19;
INDEX=floor((Sys_i_dia_mean-0.3379)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Sys_i_dia_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Sys_i_dia_mean_WEIGHT(INDEX);
end
end

%CM_mean
CM_mean_PROB=[0.7379487,0.9583875,0.94,0.8630137,0.7777778,0.8666667,0.6666667,0.25,NaN,NaN,0,NaN,NaN,NaN,NaN,0,NaN,NaN,0,NaN];
CM_mean_WEIGHT=[0.6397638,0.2522966,0.0656168,0.02395013,0.008858268,0.00492126,0.001968504,0.001312336,0,0,0.000656168,0,0,0,0,0.000328084,0,0,0.000328084,0];
DIVIDER=(0.1864-0.0002)/19;
INDEX=floor((CM_mean-0.0002)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =CM_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*CM_mean_WEIGHT(INDEX);
end
end

%SACDFMiBEM
SACDFMiBEM_PROB=[0.7246377,0.7350428,0.6629527,0.7572815,0.9282671,0.8571429,0.09090909,0.03030303,0,0,0,0,0,0,NaN,NaN,0,0,0,NaN];
SACDFMiBEM_WEIGHT=[0.0226378,0.03838583,0.1177822,0.2027559,0.4619423,0.1194226,0.0144357,0.01082677,0.004593176,0.00164042,0.001968504,0.000328084,0.000656168,0.000984252,0,0,0.000328084,0.000656168,0.000656168,0];
DIVIDER=(6.7908-0.109)/19;
INDEX=floor((SACDFMiBEM-0.109)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =SACDFMiBEM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*SACDFMiBEM_WEIGHT(INDEX);
end
end

%UACDFMiBEM
UACDFMiBEM_PROB=[0.6428571,0.6774194,0.7333333,0.7272727,0.7356322,0.6232395,0.8114856,0.9422942,0.3692308,0.25,0.06451613,0.04761905,0,0.5,0,NaN,0.6666667,0,0,NaN];
UACDFMiBEM_WEIGHT=[0.004593176,0.0101706,0.01968504,0.02526247,0.05708661,0.09317585,0.2627953,0.4662074,0.02132546,0.01574803,0.0101706,0.006889764,0.00328084,0.000656168,0.000328084,0,0.000984252,0.000328084,0.001312336,0];
DIVIDER=(4.787-0.2462)/19;
INDEX=floor((UACDFMiBEM-0.2462)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =UACDFMiBEM_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*UACDFMiBEM_WEIGHT(INDEX);
end
end

%Skew_Glob_HF
Skew_Glob_HF_PROB=[0.9321574,0.8757281,0.7584803,0.6161616,0.539823,0.527027,0.65,0.7666667,0.6333333,0.4615385,0.4,0.6,0.2,0.75,0.6666667,NaN,0.5,0.5,0,1];
Skew_Glob_HF_WEIGHT=[0.2417979,0.3379265,0.2417979,0.06496063,0.03707349,0.02427822,0.01312336,0.00984252,0.00984252,0.008530184,0.00328084,0.00164042,0.00164042,0.001312336,0.000984252,0,0.000656168,0.000656168,0.000328084,0.000328084];
DIVIDER=(41.2734-0.3297)/19;
INDEX=floor((Skew_Glob_HF-0.3297)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Glob_HF_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Glob_HF_WEIGHT(INDEX);
end
end

%Skew_Glob_UF
Skew_Glob_UF_PROB=[0.9236726,0.6331096,0.6394052,0.6337209,0.6504065,0.6666667,0.575,0.5357143,0.7857143,0.8823529,0.7777778,0.7,0.4285714,1,1,1,NaN,NaN,NaN,1];
Skew_Glob_UF_WEIGHT=[0.5931758,0.1466535,0.08825459,0.05643045,0.04035433,0.03248031,0.01312336,0.009186352,0.004593176,0.005577428,0.002952756,0.00328084,0.002296588,0.000656168,0.000328084,0.000328084,0,0,0,0.000328084];
DIVIDER=(122.6454-0.3889)/19;
INDEX=floor((Skew_Glob_UF-0.3889)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_Glob_UF_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_Glob_UF_WEIGHT(INDEX);
end
end

%Skew_S1_UF_mean
Skew_S1_UF_mean_PROB=[0.9978237,0.9531443,0.6688312,0.6145833,0.377193,0.4754098,0.5966387,0.5486111,0.5874126,0.6052632,0.6181818,0.6486486,0.7884616,0.7179487,0.7894737,0.6363636,0.5,NaN,0,1];
Skew_S1_UF_mean_WEIGHT=[0.3015092,0.2660761,0.05052494,0.03149606,0.03740158,0.04002625,0.039042,0.04724409,0.04691601,0.03740158,0.03608924,0.02427822,0.01706037,0.01279528,0.006233596,0.003608924,0.001312336,0,0.000656168,0.000328084];
DIVIDER=(6.2519-0.2799)/19;
INDEX=floor((Skew_S1_UF_mean-0.2799)/DIVIDER)+1;
if INDEX>0 & INDEX<20
MP=0;
PROB =Skew_S1_UF_mean_PROB(INDEX);
if PROB<L1
MP=-(L1-PROB);
end
if PROB>L2
MP=PROB-L2;
end
if ~isnan(MP)
GSKORE=GSKORE+MP*Skew_S1_UF_mean_WEIGHT(INDEX);
end
end

if GSKORE>-0.001540238
classifyResult=-1;
end

if GSKORE<-0.001540238
classifyResult=1;
end

%%%%%%%%% END AUTOGENERATED CODE FROM Histofind %%%%%%%%%%%%%%%%
    end
    
end

if (doPictures)
    
%  R referenceResult = Ref(2.N);
   
  FigHandle = figure('Name',recordName);
  set(FigHandle, 'Position', [100 400 2400 1000]);
  
  str = 'Abnormal';
  color = [1 0.5 0.0];
  
  if correctAnswer == -1
       str = 'Normal';
       color = [0.5 1 0];
  end
    
  
  
%    bcol = [1 1 1];
     
   %if (correctAnswer==-1)
     bcol = [0.75 1 0.95];
   %end
   
    if (correctAnswer~=classifyResult)
     bcol = [1 0.3 0.0];
   end
   
   if (classifyResult==0)
     bcol = [0.5 0.5 0.5];
   end
  
%      if (frequencyRatio<0.05)
%        bcol = [1 0.9 0.9];
%      end
   
   handle = subplot(5,1,1);
   
   
     
  
   set(handle, 'color',bcol);
    


   hold on;
   
   t1 = (1:length(PCG_resampled))./Fs;
  
   
 % y = zeros(1,length(s1positions)); %min(PCG_resampled);
   y=0;
 
   s1xPos = locExtrema(s1positions)./Fs;
   s2xPos = locExtrema(s2positions)./Fs;
   rw  = s2xPos-s1xPos;
   h = max(PCG_resampled)-min(PCG_resampled);
   
 
   if length(s1xPos)*length(s2xPos)>0
        line([s1xPos s2xPos],[y y],'color',[1,0.7,0.2],'LineWidth',40);
    end
   
   plot(t1,PCG_resampled,'k');
   
   
    if length(s1xPos)==0 | length(s2xPos)==2
       return;
    end
   
    
    
   yinvalids = max(PCG_resampled);
   invScreen = invalid'.*yinvalids;
   
   area(t1,invScreen,'facecolor',[0.7 0.7 0.7]);

%    rectangle('Position'.[s1xPos.yv.s2xPos.hv]);   

   legend('Detected S1-S2','Raw data');

   
   text (0,1,max(PCG_resampled),['FreqRatio = ' num2str(frequencyRatio)]);
   
   
   
   handle = subplot(5,1,2);
   set(handle, 'color',bcol);
     
   hold on;
   plot(t1,lowFEnvelope,'k');
%   plot(midFEnvelope,'r');
%   plot(highFEnvelope,'b');
   
   %plot(locExtrema./Fs.extrema.'o'.'color'.[0.7 0.7 0.7]);
   
   plot(locExtrema(s1positions)./Fs,extrema(s1positions),'o','MarkerFaceColor','k');
   plot(locExtrema(s2positions)./Fs,extrema(s2positions),'+','MarkerFaceColor','k');
%   plot(locExtrema(indexesToRemove)./Fs,extrema(indexesToRemove),'rs');
   
%   line([0 max(t1)],[threshold threshold]);
   legend(['Envelope LF ' num2str(lfMin) '-' num2str(lfMax) ' Hz'],'S1 ','S2');
   
   subplot(5,1,3);
    bar(fftcorrs);
    legend('Hist CorrFFT');
    
%    imagesc(spectra')
   
   subplot(5,1,4);
   hold on;
   
   xvals = 1:length(avgLShape);
   
   plot(xvals,avgLShape,'g',xvals,avgMShape.*2,'b', xvals,avgHShape.*5,'r',xvals,avgSShape.*20,'k',xvals,avgUShape.*80,'m');
   
    legend('avgLShape 15-90Hz','2 x avgMShape 55-150Hz','5 x avgHShape 100-250Hz','20 x avgSShape 200-450 Hz','80xavgUShape 400-800Hz');
   
%    plot(corels);
%    plot(sumFunc);
    
   subplot(5,1,5);
   hold on;
   
   
   plot(mdetcorrs);
   mdetratio = max(mdetcorrs)/mean(mdetcorrs);
   
   %mdetcorrs(1)=[];
   %detcorrs(end) = [];
   
   legend('Corr FFT');
   
    annotation('textbox',[0 0.95 1 0.05],'String',['Correct answer:' str],'BackgroundColor',color);
   
%   rectangle('Position',[estimatedSystoleMin,min(histogram),sr,max(histogram)],'FaceColor',[0.5 1 1]);
%   plot (histoX.histogram);
%   legend('Histogram Sys');

    color = [1 1 1];
    dim = [0.01 0.425 .1 .5];
  %  str = ['RR count =',num2str(RR_count),char(10),'RRs MEAN =',num2str(RR_mean),char(10),'RRs STD =',num2str(RR_std),char(10),'SW MEAN =',num2str(SystoleWidth_mean),char(10),'SW STD =',num2str(SystoleWidth_std)];
    
%     str = [str.char(10). 'SumFunc mean ='.num2str(mean(sumFunc)).char(10).'SumFunc STD ='.num2str(std(sumFunc)).char(10).'SumFunc P75 ='.num2str(prctile(sumFunc.75)).char(10).'SumFunc K ='.num2str(kurtosis(sumFunc))];
%     if exist('S1D')
%     str=[str,char(10),char(10),'SYM S1 COR=',num2str(mean(symmetryS1Cor)),'+-',num2str(std(symmetryS1Cor))];
%     str=[str,char(10),'SYM S1 COR MF=',num2str(mean(symmetryS1CorMF)),'+-',num2str(std(symmetryS1CorMF))];
%     str=[str,char(10),'SYM S1 COR HF=',num2str(mean(symmetryS1CorHF)),'+-',num2str(std(symmetryS1CorHF))];
%     
%     
%     str=[str,char(10),char(10),'SYM S2 COR=',num2str(mean(symmetryS2Cor)),'+-',num2str(std(symmetryS2Cor))];
%     str=[str,char(10),'SYM S2 COR MF=',num2str(mean(symmetryS2CorMF)),'+-',num2str(std(symmetryS2CorMF))];
%     str=[str,char(10),'SYM S2 COR HF=',num2str(mean(symmetryS2CorHF)),'+-',num2str(std(symmetryS2CorHF))];
%     end;
%     
%     str=[str,char(10),'SKW RAW =',num2str(skewness(PCG_resampled))];
%     str=[str,char(10),'SKW lF =',num2str(skewness(lowFEnvelope))];
%     str=[str,char(10),'SKW mF =',num2str(skewness(highFEnvelope))];
%     str=[str,char(10),'SKW hF =',num2str(skewness(midFEnvelope))];
%     
%     str=[str,char(10),char(10),'Corr Max =',num2str(max(mdetcorrs))];
%     str=[str,char(10),'Corr 75pt =',num2str(prctile(mdetcorrs,75))];
%     str=[str,char(10),'Corr Med =',num2str(median(mdetcorrs))];
%     str=[str,char(10),'Corr mean =',num2str(mean(mdetcorrs))];
%     str=[str,char(10),'Corr SKEW =',num2str(skewness(mdetcorrs))];
%     str=[str,char(10),'Corr KURT =',num2str(kurtosis(mdetcorrs))];
%     str=[str,char(10),'Corr std =',num2str(std(mdetcorrs))];
%     str=[str,char(10),'MDETratio =',num2str(mdetratio)];
    
    
    
    
    if wrongForFeatureExtraction==1
        color = [0.80 0.8 0.8];
    end
    
 %   annotation('textbox',dim,'String',str,'BackgroundColor',color);
  
    
%    text (0.1.0.);
   
   suffix ='';
   if correctAnswer==-1
       suffix = '_NORMAL';
   end
   
   
   
   hold off
   hgexport(FigHandle, [recordName,suffix,'.png'], hgexport('factorystyle'), 'Format', 'png');
   close(FigHandle);
end

if doTable
   
   logFileName = 'ver25_big3.txt';
    
   createHeaders = 0;
   fileEx = exist(logFileName, 'file');
  
   if  fileEx== 0
       createHeaders = 1;
   end
    
    
   fid = fopen(logFileName,'at');

   
   if createHeaders 
   fprintf(fid,'%s;', 'Filename');
   fprintf(fid,'%s;%s;%s;','Group','ClassifyResult','SQI');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;', 'correctAnswer','RR_count','RR_mean','RR_std','SystoleWidth_mean','SystoleWidth_std');

   fprintf(fid,'%s;%s;%s;%s;','Diastole_mean','Diastole_std','Sys_i_dia_mean','Sys_i_dia_std');

   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','AS_mean','BS_mean','CS_mean','DS_mean','ES_mean','FS_mean','GS_mean');

   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','AD_mean','BD_mean','CD_mean','DD_mean','ED_mean','FD_mean','GD_mean');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;','ADiBD','ADiCD','CDiDD','DDiFD','BDiED','EDiFD','DDiED','ADiGD');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','AM_mean','BM_mean','CM_mean','DM_mean','EM_mean','FM_mean','GM_mean');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;','AMiBM','AMiCM','BMiEM','BMiCM','DMiFM','DMiEM','AMiGM','ACDFMiBEM');
   %superF
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','SAS_mean','SBS_mean','SCS_mean','SDS_mean','SES_mean','SFS_mean','SGS_mean');

   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','SAD_mean','SBD_mean','SCD_mean','SDD_mean','SED_mean','SFD_mean','SGD_mean');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;','SADiBD','SADiCD','SCDiDD','SDDiFD','SBDiED','SEDiFD','SDDiED','SADiGD');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','SAM_mean','SBM_mean','SCM_mean','SDM_mean','SEM_mean','SFM_mean','SGM_mean');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;','SAMiBM','SAMiCM','SBMiEM','SBMiCM','SDMiFM','SDMiEM','SAMiGM','SACDFMiBEM');
   %superF
   
    %ultraF
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','UAS_mean','UBS_mean','UCS_mean','UDS_mean','UES_mean','UFS_mean','UGS_mean');

   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','UAD_mean','UBD_mean','UCD_mean','UDD_mean','UED_mean','UFD_mean','UGD_mean');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;','UADiBD','UADiCD','UCDiDD','UDDiFD','UBDiED','UEDiFD','UDDiED','UADiGD');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','UAM_mean','UBM_mean','UCM_mean','UDM_mean','UEM_mean','UFM_mean','UGM_mean');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;','UAMiBM','UAMiCM','UBMiEM','UBMiCM','UDMiFM','UDMiEM','UAMiGM','UACDFMiBEM');
   %ultraF
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;','UAMiAM','UBMiBM','UCMiCM','UDMiDM','UEMiEM','UFMiFM','UGMiGM');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;','Sym_S1_Sum_mean','Sym_S1_Sum_std','Sym_S1_std_mean', 'Sym_S1_std_std', 'Sym_S1_Corr_mean','Sym_S1_Corr_std');
   fprintf(fid,'%s;%s;%s;%s;%s;%s;','Sym_S2_Sum_mean','Sym_S2_Sum_std','Sym_S2_std_mean',' Sym_S2_std_std', 'Sym_S2_Corr_mean','Sym_S2_Corr_std');
    
   fprintf(fid,'%s;%s;%s;%s;','Sym_S1_Sum_mean_i_std','Sym_S2_Sum_mean_i_std','Sym_S1_Corr_mean_i_std','Sym_S2_Corr_mean_i_std');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;','Sym_S1_Corr_MF_mean','Sym_S1_Corr_MF_std','Sym_S2_Corr_MF_mean','Sym_S2_Corr_MF_std','Sym_S1_Corr_HF_mean','Sym_S1_Corr_HF_std','Sym_S2_Corr_HF_mean','Sym_S2_Corr_HF_std','Sym_S1_Corr_SF_mean','Sym_S1_Corr_SF_std','Sym_S2_Corr_SF_mean','Sym_S2_Corr_SF_std','Sym_S1_Corr_UF_mean','Sym_S1_Corr_UF_std','Sym_S2_Corr_UF_mean','Sym_S2_Corr_UF_std');

   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;','Sym_S1_Corr_MF_mean_i_std','Sym_S2_Corr_MF_mean_i_std','Sym_S1_Corr_HF_mean_i_std','Sym_S2_Corr_HF_mean_i_std','Sym_S1_Corr_SF_mean_i_std','Sym_S2_Corr_SF_mean_i_std','Sym_S1_Corr_UF_mean_i_std','Sym_S2_Corr_UF_mean_i_std');
   
   fprintf(fid,'%s;%s;%s;','Sym_S1_Corr_MF_std_i_Sym_S2_Corr_MF_mean','Sym_S1_Corr_mean_i_Sym_S2_Corr_mean','Sym_S1_Corr_std_i_Sym_S2_Corr_std');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;','Skew_Glob_Raw','Skew_Glob_LF','Skew_Glob_MF','Skew_Glob_HF','Skew_Glob_SF','Skew_Glob_UF');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;','Skew_Loc_Raw_mean','Skew_Loc_Raw_std','Skew_Loc_LF_mean','Skew_Loc_LF_std','Skew_Loc_MF_mean','Skew_Loc_MF_std','Skew_Loc_HF_mean','Skew_Loc_HF_std','Skew_Loc_SF_mean','Skew_Loc_SF_std','Skew_Loc_UF_mean','Skew_Loc_UF_std');
   
   fprintf(fid,'%s;%s;%s;%s;','Skew_S1_LF_mean','Skew_S1_LF_std','Skew_S2_LF_mean','Skew_S2_LF_std');
   fprintf(fid,'%s;%s;%s;%s;','Skew_S1_MF_mean','Skew_S1_MF_std','Skew_S2_MF_mean','Skew_S2_MF_std');
   fprintf(fid,'%s;%s;%s;%s;','Skew_S1_HF_mean','Skew_S1_HF_std','Skew_S2_HF_mean','Skew_S2_HF_std');
   fprintf(fid,'%s;%s;%s;%s;','Skew_S1_SF_mean','Skew_S1_SF_std','Skew_S2_SF_mean','Skew_S2_SF_std');
   fprintf(fid,'%s;%s;%s;%s;','Skew_S1_UF_mean','Skew_S1_UF_std','Skew_S2_UF_mean','Skew_S2_UF_std');
   
   fprintf(fid,'%s;%s;%s;%s;','jhParam_mean','jhParam_std','jhParam_skew','jhParam_kurt');
  
   fprintf(fid,'%s;%s;%s;%s;%s;','XT_L','XT_M','XT_H','XT_S','XT_U');
   
   fprintf(fid,'%s;%s;%s;%s;%s;%s;','Corr_FFT_mean','Corr_FFT_std','Corr_FFT_skw','Corr_FFT_kurt','Corr_FFT_median','Corr_FFT_mean_i_std');

    fprintf(fid,'%s;%s;%s;%s;%s;%s;','DCODE','limitCheck','frequencyRatio','invalidPart','globalHFratio','DIAG');

   
   fprintf(fid,'%s\n', ' ');
   end
   
if true %exist('RR_count')   
   
   str = [recordName,';'];
   fprintf(fid,'%s;', str);
   fprintf(fid,'%1.0f;%1.0f;%1.0f;',class,classifyResult,quality);
   fprintf(fid,'%1.0f;', correctAnswer);

   
  

if (minHistoVal~=maxHistoVal) & (length(s1positions)>2) 
   
   fprintf(fid,'%1.0f;%1.10f;%1.10f;%1.10f;%1.10f;',RR_count,RR_mean,RR_std,SystoleWidth_mean,SystoleWidth_std)
    
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',Diastole_mean,Diastole_std,Sys_i_dia_mean,Sys_i_dia_std);
    
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',AS_mean,BS_mean,CS_mean,DS_mean,ES_mean,FS_mean,GS_mean);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',AD_mean,BD_mean,CD_mean,DD_mean,ED_mean,FD_mean,GD_mean);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',ADiBD,ADiCD,CDiDD,DDiFD,BDiED,EDiFD,DDiED,ADiGD);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',AM_mean,BM_mean,CM_mean,DM_mean,EM_mean,FM_mean,GM_mean);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',AMiBM,AMiCM,BMiEM,BMiCM,DMiFM,DMiEM,AMiGM,ACDFMiBEM);
   
   %super F
     fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',SAS_mean,SBS_mean,SCS_mean,SDS_mean,SES_mean,SFS_mean,SGS_mean);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',SAD_mean,SBD_mean,SCD_mean,SDD_mean,SED_mean,SFD_mean,SGD_mean);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',SADiBD,SADiCD,SCDiDD,SDDiFD,SBDiED,SEDiFD,SDDiED,SADiGD);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',SAM_mean,SBM_mean,SCM_mean,SDM_mean,SEM_mean,SFM_mean,SGM_mean);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',SAMiBM,SAMiCM,SBMiEM,SBMiCM,SDMiFM,SDMiEM,SAMiGM,SACDFMiBEM);
   %%%
   
     %ultra F
     fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',UAS_mean,UBS_mean,UCS_mean,UDS_mean,UES_mean,UFS_mean,UGS_mean);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',UAD_mean,UBD_mean,UCD_mean,UDD_mean,UED_mean,UFD_mean,UGD_mean);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',UADiBD,UADiCD,UCDiDD,UDDiFD,UBDiED,UEDiFD,UDDiED,UADiGD);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',UAM_mean,UBM_mean,UCM_mean,UDM_mean,UEM_mean,UFM_mean,UGM_mean);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',UAMiBM,UAMiCM,UBMiEM,UBMiCM,UDMiFM,UDMiEM,UAMiGM,UACDFMiBEM);
   %%%
   
     fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',UAMiAM,UBMiBM,UCMiCM,UDMiDM,UEMiEM,UFMiFM,UGMiGM);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',Sym_S1_Sum_mean,Sym_S1_Sum_std,Sym_S1_std_mean, Sym_S1_std_std, Sym_S1_Corr_mean,Sym_S1_Corr_std);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',Sym_S2_Sum_mean,Sym_S2_Sum_std,Sym_S2_std_mean, Sym_S2_std_std, Sym_S2_Corr_mean,Sym_S2_Corr_std);
    
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',Sym_S1_Sum_mean_i_std,Sym_S2_Sum_mean_i_std,Sym_S1_Corr_mean_i_std,Sym_S2_Corr_mean_i_std);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',Sym_S1_Corr_MF_mean,Sym_S1_Corr_MF_std,Sym_S2_Corr_MF_mean,Sym_S2_Corr_MF_std,Sym_S1_Corr_HF_mean,Sym_S1_Corr_HF_std,Sym_S2_Corr_HF_mean,Sym_S2_Corr_HF_std,Sym_S1_Corr_SF_mean,Sym_S1_Corr_SF_std,Sym_S2_Corr_SF_mean,Sym_S2_Corr_SF_std,Sym_S1_Corr_UF_mean,Sym_S1_Corr_UF_std,Sym_S2_Corr_UF_mean,Sym_S2_Corr_UF_std);

   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',Sym_S1_Corr_MF_mean_i_std,Sym_S2_Corr_MF_mean_i_std,Sym_S1_Corr_HF_mean_i_std,Sym_S2_Corr_HF_mean_i_std,Sym_S1_Corr_SF_mean_i_std,Sym_S2_Corr_SF_mean_i_std,Sym_S1_Corr_UF_mean_i_std,Sym_S2_Corr_UF_mean_i_std);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;',Sym_S1_Corr_MF_std_i_Sym_S2_Corr_MF_mean,Sym_S1_Corr_mean_i_Sym_S2_Corr_mean,Sym_S1_Corr_std_i_Sym_S2_Corr_std);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',Skew_Glob_Raw,Skew_Glob_LF,Skew_Glob_MF,Skew_Glob_HF,Skew_Glob_SF,Skew_Glob_UF);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',Skew_Loc_Raw_mean,Skew_Loc_Raw_std,Skew_Loc_LF_mean,Skew_Loc_LF_std,Skew_Loc_MF_mean,Skew_Loc_MF_std,Skew_Loc_HF_mean,Skew_Loc_HF_std,Skew_Loc_SF_mean,Skew_Loc_SF_std,Skew_Loc_UF_mean,Skew_Loc_UF_std);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',Skew_S1_LF_mean,Skew_S1_LF_std,Skew_S2_LF_mean,Skew_S2_LF_std);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',Skew_S1_MF_mean,Skew_S1_MF_std,Skew_S2_MF_mean,Skew_S2_MF_std);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',Skew_S1_HF_mean,Skew_S1_HF_std,Skew_S2_HF_mean,Skew_S2_HF_std);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',Skew_S1_SF_mean,Skew_S1_SF_std,Skew_S2_SF_mean,Skew_S2_SF_std);
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',Skew_S1_UF_mean,Skew_S1_UF_std,Skew_S2_UF_mean,Skew_S2_UF_std);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',XT_L,XT_M,XT_H,XT_S,XT_U);
   
   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;',jhParam_mean,jhParam_std,jhParam_skew,jhParam_kurt);
   
end

   fprintf(fid,'%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;%1.10f;',Corr_FFT_mean,Corr_FFT_std,Corr_FFT_skw,Corr_FFT_kurt,Corr_FFT_median,Corr_FFT_mean_i_std);

   fprintf(fid,'%1.0f;%1f.0;%1.12f;%1.8f;%1.12f;%s;',diagnoseCode,limitCheck,frequencyRatio,invalidPart,globalHFratio,'blah');
   
   fprintf(fid,'%s\n', ' ');
end    

   fclose(fid);
end


