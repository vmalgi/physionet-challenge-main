function [s1positions s2positions extrema locExtrema minHistoVal maxHistoVal histogram avgLShape avgMShape avgHShape avgSShape avgUShape corrs corrHist] = getS1S2( PCG_resampled,Fs,lowFEnvelope,midFEnvelope,highFEnvelope,superFEnvelope,ultraFEnvelope, invalid )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

dl = length(PCG_resampled);
f = 10; %z 20
[b,a] = butter(3, 2*f/Fs, 'low'); 
lowFE_filtered = filtfilt(b,a,lowFEnvelope); 
midFE_filtered = filtfilt(b,a,midFEnvelope); 
highFE_filtered = filtfilt(b,a,highFEnvelope); 
superFE_filtered = filtfilt(b,a,superFEnvelope); 
ultraFE_filtered = filtfilt(b,a,ultraFEnvelope); 

%-----------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------Estimated systole length----------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------------------------------


threshold = prctile(lowFE_filtered(find(invalid==0)),74);

%if (dl<10000) threshold=prctile(lowFE_filtered(find(invalid==0)),60);
%end

 
[extrema,locExtrema] = findpeaks(lowFE_filtered.*(1-invalid),'MinPeakHeight',threshold,'MinPeakDistance',constants.minSysDuration);

distances = diff(locExtrema);


distances (distances>constants.maxSysDuration)=[];

histoBins = 10; %nemìlo by to být proporcionální k poètu samplù? TOSEARCH

if dl<10000
    histoBins = 5;
end

maxHistoVal = prctile(distances,100);
minHistoVal = prctile(distances,0);

distances(distances<minHistoVal)=[];
distances(distances>maxHistoVal)=[];

histogram = hist(distances,histoBins);

histoStep = (maxHistoVal-minHistoVal)/(histoBins-1);
histoX = minHistoVal:histoStep:maxHistoVal;

%vyhledání maxim. Pokud to lepší významnì dominuje, pak není co øešit.
%Jinak se vybere to nižší

[histoMaxValue histoMaxPos] = findpeaks(histogram, 'SortStr','descend');

%estimatedSystoleUncertainity = 0;

[globHmaxVal globHmaxPos]  = max(histogram);



if globHmaxVal>0.7*max(histoMaxValue) % nejvìtší hodnota je okrajová
    
    if globHmaxPos==1
        histoMaxValue = [globHmaxVal histoMaxValue];
        histoMaxPos = [1 histoMaxPos];
    end
    
    if globHmaxPos==histoBins
        histoMaxValue = [globHmaxVal histoMaxValue];
        histoMaxPos = [histoBins histoMaxPos];
    end

    
end



if length(histoMaxValue)==0 %pokud je extrém na kraji
    hindex = 1;
    
    if histogram(histoBins)*0.7>histogram(1)
        hindex = histoBins;
    end
    
end

if length(histoMaxValue)==1 
        hindex = histoMaxPos(1);
end

if length(histoMaxValue)>1
    hindex = histoMaxPos(1);
    
    %-------------------------------  pišvejc 0.75------------------
    if histoMaxValue(2)>0.7*histoMaxValue(1) & histoMaxPos(2)<histoMaxPos(1) 
%        valueHMax = histoMaxValue(2);
        hindex = histoMaxPos(2);
    end
end


if (length(histoX)>0)
estimatedSystoleDuration = histoX(hindex);
else
estimatedSystoleDuration = minHistoVal; %----------nouzovka. Možná vyhodit
end



%-----------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------RR estimation--------------------------------------------------------------- 
%-----------------------------------------------------------------------------------------------------------------------------------------

autoCorr= zeros(dl,1);
step = 50; %orig 20
for i=2:step:length(autoCorr)
    weight = min([(length(autoCorr)-i)/500 i/2000 1]);
    block1 = lowFEnvelope(1:i);
    block2 = lowFEnvelope(length(autoCorr)-i+1:length(autoCorr));

    crc = corrcoef(block1,block2);
    
    autoCorr(i:i+step) = crc(1,2)*weight;
end


%f = 5; %z 20
%[b,a] = butter(3, 2*f/Fs, 'low'); 
autoCorr = fliplr(lowFE_filtered); 

%odhad peaks number
expPeaksNum = fix(dl/1000);
%ac_threshold = max(autoCorr)/4;
ac_threshold = prctile(autoCorr,80);
[rhp,rhpLocs] = findpeaks(autoCorr,'MinPeakHeight',ac_threshold, 'MinPeakDistance',constants.minRRDuration,'SortStr','descend','NPeaks',expPeaksNum);

sortedLocs = sortrows(rhpLocs);

difLocs = diff(sortedLocs);
difLocs(difLocs>constants.maxRRDuration)=[];

estimatedRRinterval = median(difLocs);  
estimatedRRintMean = mean(difLocs);
estimatedRRintSTD = std(difLocs);
estimatedRRcorr = median(rhp);


if length(difLocs)<3
    estimatedRRinterval = 2.8*estimatedSystoleDuration;
end

rrSystoleRatio = estimatedRRinterval/estimatedSystoleDuration;

%-----------------------------------------------------------------------------------------------------------------------------------------
%-------------------------------------------------------------RR/SYS korekce-------------------------------------------------------------- 
%-----------------------------------------------------------------------------------------------------------------------------------------

if rrSystoleRatio<constants.rrSysRatioLowLimit  || rrSystoleRatio>constants.rrSysRatioHighLimit

corrExpectedSysDuration = estimatedRRinterval/2.5;
systoleSolutions = abs(histoX-corrExpectedSysDuration);

[minimums minimLocs] = min(systoleSolutions);

if length(minimLocs)>0
if histogram(minimLocs(1))>expPeaksNum/3 
    estimatedSystoleDuration = histoX(minimLocs(1));
else
    if histogram(hindex) > expPeaksNum/2
    estimatedRRinterval = estimatedRRinterval+estimatedSystoleDuration;
    end
end


end
end

rrSystoleRatio = estimatedRRinterval/estimatedSystoleDuration;

if rrSystoleRatio<constants.rrSysRatioLowLimit
    
    
    if estimatedRRinterval/estimatedRRintSTD<constants.reliableRRsSTD %RR je asi považován mylnì za DS
        estimatedRRinterval = estimatedRRinterval + estimatedSystoleDuration;
    else
        estimatedSystoleDuration = estimatedRRinterval - estimatedSystoleDuration;
    end
   
    
end

if rrSystoleRatio>constants.rrSysRatioHighLimit
    if estimatedRRinterval>1000
        estimatedRRinterval = estimatedRRinterval/2; %trochu fake
    end
    
end

rrSystoleRatio = estimatedRRinterval/estimatedSystoleDuration;

estimatedSystoleMin = max(estimatedSystoleDuration-2.5*histoStep, constants.minSysDuration);
estimatedSystoleMax = min(estimatedSystoleDuration+3.5*histoStep, constants.maxSysDuration); %bylo 1.5

%estimatedSystoleMin = estimatedSystoleDuration-2.5*histoStep;
%estimatedSystoleMax = estimatedSystoleDuration+2.5*histoStep; %bylo 1.5


%----------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------S1/S2-search-------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------

s1positions = [];
s2positions = [];

%threshold = prctile(lowFE_filtered(find(invalid==0)),50);
[extrema,locExtrema] = findpeaks(lowFE_filtered.*(1-invalid),'MinPeakHeight',threshold,'MinPeakDistance',constants.minSysDuration);



for s1index = 1:length(extrema)-1
    %kontrola zda se od S1 nachází nìkde vpravo peak v rozumné vzdálenosti
    
    if length(s2positions)>=1 & s2positions(end)==s1index
        continue;
    end
    
    s1pos = locExtrema(s1index);
    s1amp = lowFE_filtered(s1pos);
    
    beforeS1start = fix(max([s1pos-estimatedSystoleDuration 1]));
    beforeS1end = fix(max([s1pos-100 1]));
    
    if beforeS1start ==1
        continue;
    end
    
    maximumFromS1ToLeft = max(lowFE_filtered(beforeS1start:beforeS1end));
    
    %S1 musí být silnìjší než oblast vlevo
    if maximumFromS1ToLeft>s1amp
        continue;
    end
    
    acceptS2=0;

    
    %najít potenciální S2 podle délky
    acceptableS2indexes = [];
    for s2index = s1index+1:length(extrema)
        s2pos = locExtrema(s2index);
        distance = s2pos-s1pos;
        if distance<estimatedSystoleMin
            continue;
        end
        
        if distance>estimatedSystoleMax
            break;
        end
        
        acceptableS2indexes = [acceptableS2indexes s2index];
    end
    
    %kontrola na nenalezení S2 peaku
    if length(acceptableS2indexes)==0
        continue;
    end
    
    % vybrat nejpravdìpodobnìjší S2 podle nejvìtší amplitudy. To ale nemusí
    % být pravda. Mìlo by se to zmenšovat s vìtším rozdílem od
    % pøedpokládané délky systoly
    %[s2amp s2indexInAccepted] = max(extrema(acceptableS2indexes));
     
    
    [s2amp s2indexInAccepted] = max(extrema(acceptableS2indexes));
    
    s2index = acceptableS2indexes(s2indexInAccepted);
    
    s2pos = locExtrema(s2index);
    
    s1s2distance = s2pos-s1pos;
    
    %s1 amp musí být vìtší než 0.5 s2amp 
    if s1amp<0.2*s2amp %pùvodnì 0.5
        continue;
    end
    
    
    %S2 musí být silnìjší než oblast vpravo
    
    afterS2start =fix(min([s2pos+100 dl]));
    afterS2end = fix(min([s2pos+s1s2distance dl]));
    
  %vynásobení oknem. Možná bude dìlat blbosti!!!
    
    
    smdata = lowFE_filtered(afterS2start:afterS2end-1);
    
    if estimatedSystoleDuration<2500 %test
        window = hamming(afterS2end-afterS2start);    
        smdata = smdata.*window;
    end
    
 %    maximumFromS2ToRight = max(smdata.*window);
     maximumFromS2ToRight = max(smdata);
    
    if s2index<length(locExtrema)
        
        locExtNextIndex = s2index+1;
        locExtNextPos = locExtrema(locExtNextIndex);
        locExtNextAmp = lowFE_filtered(locExtNextPos);
        
        distNext = locExtNextPos-s2pos;
        
        if distNext<s1s2distance & distNext>estimatedSystoleMin & locExtNextAmp>s1amp*0.6
            continue;
        end
        
    end
    
    
    
    if maximumFromS2ToRight>s2amp
        continue;
    end
    
   %kontrola na to, zda nesdílím nový s2 se starým s2
   if length(s2positions)>0 & s2index==s2positions(end)
       
   %    preS1amp = extrema(s1positions(end));
   %    newS1amp = extrema(s1index);
   %    S2amp = extrema(s2index);
       
   %    if newS1amp>preS1amp
           s1positions(end)=[];
           s2positions(end)=[];
    %   else
     %      continue;
     %  end
       
   end
    
    s1positions = [s1positions s1index];
    s2positions = [s2positions s2index];
   
    
end

if length(s1positions)*length(s2positions)==0
    avgLShape = zeros(100,1);
    avgMShape = zeros(100,1);
    avgHShape = zeros(100,1);
    avgSShape = zeros(100,1);
    avgUShape = zeros(100,1);
    corrs = zeros(100,1);
    corrHist = zeros(100,1);
else
%----------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------Systole AVG shape---------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------

systoles = locExtrema(s2positions)-locExtrema(s1positions);

meanSysDuration = mean(systoles); %nemìl by být medián????

avgRad = fix(meanSysDuration*2);
avgLength = avgRad*2;

avgLShape = zeros(avgLength,1);
avgMShape = zeros(avgLength,1);
avgHShape = zeros(avgLength,1);
avgSShape = zeros(avgLength,1);
avgUShape = zeros(avgLength,1);

globalLFmax = max(lowFE_filtered);

for pos=1:length(s1positions)
    
    s1index = s1positions(pos);
    s2index = s2positions(pos);
    
    s1pos = locExtrema(s1index);
    s2pos = locExtrema(s2index);
    
    center = (s2pos+s1pos)/2;
    left = fix(center-avgRad);
    right =fix(center+avgRad-1);
    
    if left<1 continue;
    end
    if right>dl break;
    end
    
    originalData = lowFE_filtered(left:right);
    originalDataM = midFE_filtered(left:right);
    originalDataH = highFE_filtered(left:right);
    originalDataS = superFE_filtered(left:right);
    originalDataU = ultraFE_filtered(left:right);
    
    if max(originalData)==globalLFmax
        continue;
    end
    
    avgLShape=avgLShape+originalData; 
    avgMShape = avgMShape+originalDataM;
    avgHShape = avgHShape+originalDataH;
    avgSShape = avgSShape+originalDataS;
    avgUShape = avgUShape+originalDataU;
end


%----------------------------------------------------------------------------------------------------------------------------------------------
%---------------------------------------------------------------------S1/S2-corelation add-on---------------------------------------------------
%----------------------------------------------------------------------------------------------------------------------------------------------


distances = diff(locExtrema(s1positions));

holeLimit = 1.8*estimatedRRinterval;

if max(distances)>holeLimit | locExtrema(s1positions(1))>holeLimit | dl-locExtrema(s2positions(end))>holeLimit %provádím, pouze pokud je tam díra vìtší jak 2*estRR
    
corrStart = avgRad+1;
corrEnd = dl-avgRad-1;
    
corrCurve = zeros(dl,1);    
    
% sestavení korelaèní funkce
for cpos = corrStart:2:corrEnd
    
    
    closeToAny=0;
    for pos=s2positions
        tp=locExtrema(pos);
        dist = abs(tp-cpos); 
        if dist<estimatedRRinterval*0.7
            closeToAny = 1;
            break;
        end
    end
    
    if closeToAny
        continue;
    end
    
    
    windowStart = fix(cpos-avgRad);
    windowEnd = windowStart+length(avgLShape)-1;
    
    lfdata = lowFE_filtered(windowStart:windowEnd);
    
     crc = corrcoef(lfdata,avgLShape);
        
     corrCurve(cpos)=crc(1,2);
end
    

% dosazení chybìjících peakù

threshCorr = 0.8;
if max(corrCurve>=threshCorr)

[peakAmps,locExtraPeaks] = findpeaks(corrCurve,'MinPeakHeight',threshCorr,'MinPeakDistance',fix(estimatedSystoleMax));    

[mv leftOffsetS1] = max(avgLShape(fix(avgRad/2):avgRad));
[mv rightOffsetS2] =max(avgLShape(avgRad:fix(avgRad+avgRad/2)));

leftOffsetS1 =  leftOffsetS1-fix(avgRad/2);

locExtrema = locExtrema';
extrema = extrema';

for extraCenter = 1:length(locExtraPeaks)
    centerPos = locExtraPeaks(extraCenter);
    
    s1pos = centerPos+leftOffsetS1;
    s2pos = centerPos+rightOffsetS2;
    
  locExtrema = [locExtrema s1pos];
  extrema = [extrema 0];
  s1positions = [s1positions length(locExtrema)];
  
  locExtrema = [locExtrema s2pos];
  extrema = [extrema 0];
  s2positions = [s2positions length(locExtrema)];
end


locExtrema = locExtrema';
extrema = extrema';
end

% seøazení
if true

finalLocExtrema = [];
finalS1positions = [];
finalS2positions = [];

finalPoints = [];

lp = length(s1positions);

for s1ind = 1:length(s1positions)
    finalPoints(s1ind,1)=locExtrema(s1positions(s1ind));
    finalPoints(s1ind,2)=locExtrema(s2positions(s1ind));
end

finalPoints = sort(finalPoints);

locExtrema = [];
s1positions=[];
s2positions = [];
for s1ind = 1:lp
    
    locExtrema(end+1) = finalPoints(s1ind,1);
    s1positions(s1ind)=length(locExtrema);
    
    locExtrema(end+1) = finalPoints(s1ind,2);
    s2positions(s1ind)=length(locExtrema);
end

locExtrema = locExtrema';
extrema = lowFE_filtered(locExtrema);
s1positions = s1positions';
s2positions = s2positions';

end


end

end


%--------------------------------------MURMURTEST---------
%---------------------------------------------------------

winsize = 128;

corrs = zeros(dl,1);

if true

off = 1;
    
lastSpectra = zeros(winsize/2-off,1);


registeredDelays =  zeros(winsize/2,1);
lastExcitedSamplePos =  zeros(winsize/2,1);
numColumns = zeros(winsize/2,1);

window = hamming(winsize+1);

step = winsize/16;

thresholds = zeros(winsize/2,fix((dl-winsize)/step)+1);

for pos=1:step:dl-winsize
    testedWindow=PCG_resampled(pos:pos+winsize).*window;
    
    
    
    fftimage = fft(testedWindow,winsize);
    
    
    
    spectra = abs(fftimage(1+off:winsize/2));
    
    %spectra = real(fft(testedWindow,winsize));
    
     crc = corrcoef(spectra,lastSpectra);
        
     corrs(pos:pos+step)=crc(1,2)*sum(log(spectra));
    
    lastSpectra = spectra;
end





corrs(1)=[];
corrs(dl-winsize:end)=[];


corrs(isnan(corrs)) = [] ;

 f = 0.5; %z 20
 [b,a] = butter(3, 2*f/Fs, 'high'); 
 corrs = filtfilt(b,a,corrs); 


corrHist = hist(corrs,fix(winsize/2-off));

% 
% fftim2 = abs(fft(corrs,length(corrs)));
% corFFTim = log(fftim2(fix(1:length(fftim2)/2)));
% 


end




%toto bylo použité døív. Nebylo to úplnì blbé, ale dìlalo to chyby.
%distances = diff(locExtrema);
%s1positions = find((distances>estimatedSystoleMin) & (distances<estimatedSystoleMax));
%s2positions = s1positions+1;

end

