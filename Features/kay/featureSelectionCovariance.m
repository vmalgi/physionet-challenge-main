%This script tries to take a normalized features matrix and reduces its
%size by removing highly correlated variables

clear all
sampling = 'custom';
featuresString = 'holidayTest12AllNormalizedtTest585';
covarianceThreshold = 0.9;
removedFeatures = [];

%% %%%%%%%%%%%%%%%%%%%%%%  Load in features matrix for training data %%%%%%%%%%%%%%%%%%%%%


NoValLoad = strcat('Features',filesep,'NOVAL',sampling,featuresString,'.mat');
load(NoValLoad)
if ~exist('featuresMatrixTrain')
    if size(featuresMatrix,1) > size(featuresMatrix,2)
        featuresMatrixTrain = featuresMatrix';
    else
        featuresMatrixTrain = featuresMatrix;
    end
end

%% %%%%%%%%%%%%% Get covariance matrix and work out which features to remove %%%%%%%%%%% 
covarianceMatrix = cov(featuresMatrixTrain');


for i = 1:1000
    for i = 1:size(covarianceMatrix,1);
        covarianceMatrix(i,i) = 0;
    end
    [maxCovariance,indexRow] = max(abs(covarianceMatrix));
    [sortedMaxCovariances,iSort] = sort(maxCovariance,'descend');
    indexRow = indexRow(iSort);
    
    if sortedMaxCovariances(1) > covarianceThreshold
        
        %pick the one with the highest value on a tTest
        featuresMatrixNormal = featuresMatrixTrain(:,find(outcomes==-1));
        featuresMatrixAbnormal = featuresMatrixTrain(:,find(outcomes==1));
        
        
        featuresMeanNormal = mean(featuresMatrixNormal,2);
        featuresVarNormal = var(featuresMatrixNormal,0,2);
        
        featuresMeanAbnormal = mean(featuresMatrixAbnormal,2);
        featuresVarAbnormal = var(featuresMatrixAbnormal,0,2);
        
        nNormal = size(featuresMatrixNormal,2);
        nAbnormal = size(featuresMatrixAbnormal,2);
        
        sx1Minusx2 = sqrt((featuresVarNormal./nNormal) + ((featuresVarAbnormal./nAbnormal)));
        t = abs((featuresMeanNormal - featuresMeanAbnormal)./(sx1Minusx2));
        
        
        if t(indexRow(1)) < t(indexRow(2))
            covarianceMatrix(indexRow(1),:) = zeros(1,size(featuresMatrixTrain,1));
            covarianceMatrix(:,indexRow(1)) = zeros(size(featuresMatrixTrain,1),1);
            removedFeatures = [removedFeatures;indexRow(1)]
        else
            covarianceMatrix(indexRow(2),:) = zeros(1,size(featuresMatrixTrain,1));
            covarianceMatrix(:,indexRow(2)) = zeros(size(featuresMatrixTrain,1),1);
            removedFeatures = [removedFeatures;indexRow(2)]
        end
    end
end



First = 226;
Second = First + 275;
Third = Second + 15;
Fourth = Third + 0;
Fifth = Fourth + 21;
Sixth = Fifth + 22;
Seventh = Sixth + 17;
% Eighth = Seventh + 15;

RemovedCWT = sum(removedFeatures <= First)
RemovedCepstrum = sum(removedFeatures > First & removedFeatures <= Second)
RemovedInterBeat = sum(removedFeatures > Second & removedFeatures <= Third)
RemovedZCP = sum(removedFeatures > Third & removedFeatures <= Fourth)
RemovedSimplicity = sum(removedFeatures > Fourth & removedFeatures <= Fifth)
RemovedApEn = sum(removedFeatures > Fifth & removedFeatures <= Sixth)
RemovedSpecS = sum(removedFeatures > Sixth & removedFeatures <= Seventh)
RemovedSSK = sum(removedFeatures>Seventh)


featuresMatrixTrain(removedFeatures,:) = [];

NoValSave = strcat('Features',filesep,'NOVAL',sampling,featuresString,'CovarianceGT',num2str(covarianceThreshold*100),'Removed','.mat');
save(NoValSave,'featuresMatrixTrain','outcomes','signalQuality')


%% %%%%%%%%%%%%%%%%Load in validation data and remove the same feature %%%%%%%%%%%%%%%%%%%%%%%%
clearvars outcomes
clearvars signalQuality
ValLoad = strcat('Features',filesep,'Val',sampling,featuresString,'.mat');
load(ValLoad)
featuresMatrixVal(removedFeatures,:) = [];
ValSave = strcat('Features',filesep,'Val',sampling,featuresString,'CovarianceGT',num2str(covarianceThreshold*100),'Removed','.mat');
% ValSave = strcat('Features',filesep,'Val',sampling,'holidayTest1','CovarianceGT',num2str(covarianceThreshold*100),'Removed','.mat');
save(ValSave,'featuresMatrixVal','outcomes','signalQuality')




