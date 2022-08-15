%Ensemble Network
clear all
close all

profile on

featuresString = '200features50300CWT3737';
NoValLoad = strcat('Features',filesep,'NOVAL',featuresString,'.mat');
load(NoValLoad)

featuresMatrixTrain = featuresMatrix';
outcomesTrain = outcomes';

for i = 1:length(outcomesTrain)
    if outcomesTrain(i) == -1
        outcomesTrain(i) = 0;
    end
end


ValLoad = strcat('Features',filesep,'Val',featuresString,'.mat');
load(ValLoad)

featuresMatrixVal = featuresMatrix';
outcomesVal = outcomes';



% hiddenNodesVector = [70 100];
% epochsVector = [50 90];
% miniBatchSizeVector = [8];
% learningRateVector = [0.3];
% lmdbaVector = [0.2 0.3];
% momentumVector = [0.2 0.4];
% nRepeats = 7;
%networks = {'7811-Apr-2016 18:19:38220featuresEnvelope50300CWT3737','7811-Apr-2016 20:38:40220featuresEnvelope50300CWT3737','7811-Apr-2016 20:44:00220featuresEnvelope50300CWT3737','7812-Apr-2016 03:34:40220featuresEnvelope50300CWT3737','7712-Apr-2016 01:26:15220featuresEnvelope50300CWT3737','7712-Apr-2016 00:48:46220featuresEnvelope50300CWT3737','7712-Apr-2016 00:42:41220featuresEnvelope50300CWT3737','7712-Apr-2016 00:38:54220featuresEnvelope50300CWT3737','7711-Apr-2016 22:14:53220featuresEnvelope50300CWT3737','7711-Apr-2016 21:18:59220featuresEnvelope50300CWT3737','7711-Apr-2016 20:33:19220featuresEnvelope50300CWT3737','7711-Apr-2016 19:00:52220featuresEnvelope50300CWT3737','7711-Apr-2016 18:08:44220featuresEnvelope50300CWT3737','7711-Apr-2016 18:07:03220featuresEnvelope50300CWT3737'};
%networks = {'7713-Apr-2016 18:19:01200features50300CWT3737','7713-Apr-2016 20:58:34200features50300CWT3737','7713-Apr-2016 21:36:01200features50300CWT3737','7714-Apr-2016 09:15:45200features50300CWT3737','7714-Apr-2016 09:30:48200features50300CWT3737','7714-Apr-2016 09:32:42200features50300CWT3737','7714-Apr-2016 09:44:03200features50300CWT3737','7714-Apr-2016 09:56:34200features50300CWT3737','7714-Apr-2016 10:18:52200features50300CWT3737','7714-Apr-2016 11:50:17200features50300CWT3737','7714-Apr-2016 12:15:07200features50300CWT3737','7714-Apr-2016 15:02:48200features50300CWT3737','7813-Apr-2016 18:05:46200features50300CWT3737','7813-Apr-2016 18:49:51200features50300CWT3737','7814-Apr-2016 09:35:16200features50300CWT3737','7814-Apr-2016 10:00:47200features50300CWT3737','7814-Apr-2016 12:56:46200features50300CWT3737','7814-Apr-2016 13:06:28200features50300CWT3737','7814-Apr-2016 13:30:30200features50300CWT3737','7814-Apr-2016 13:42:56200features50300CWT3737','7913-Apr-2016 23:19:19200features50300CWT3737','8014-Apr-2016 14:03:42200features50300CWT3737'};
%networks = {'7813-Apr-2016 18:05:46200features50300CWT3737','7813-Apr-2016 18:49:51200features50300CWT3737','7814-Apr-2016 09:35:16200features50300CWT3737','7814-Apr-2016 10:00:47200features50300CWT3737','7814-Apr-2016 12:56:46200features50300CWT3737','7814-Apr-2016 13:06:28200features50300CWT3737','7814-Apr-2016 13:30:30200features50300CWT3737','7814-Apr-2016 13:42:56200features50300CWT3737','7913-Apr-2016 23:19:19200features50300CWT3737','8014-Apr-2016 14:03:42200features50300CWT3737'};
load('GoodNetworks/bestNetworksEnsemble')
networks = networksString;
size(networks)
nNetworksInEnsemble = [7];
nPicksPerNEnsemble = 1;
percentagesMatrix = [];
networksUsed = cell(length(nNetworksInEnsemble).*nPicksPerNEnsemble,1);
N = 0;
for jj = 1:length(nNetworksInEnsemble)
    for kk = 1:nPicksPerNEnsemble
        networksChosen = datasample(networks,nNetworksInEnsemble(jj),'Replace',false);
        eachNetworkGuess = zeros(length(outcomesVal),length(networks));
        finalGuess = zeros(length(outcomesVal),1);
        N = N + 1
        for oo = 1:length(networksChosen)
            clearvars net
            load(strcat('GoodNetworks',filesep,networksChosen{oo},'.mat'))
            
            for i = 1:length(outcomesVal)
                predictor = feed_forward(net,featuresMatrixVal(:,i));
                thr  = 0.5; % classification threshold, thr>0.5 for abnormal recordings.
                if predictor > thr
                    classifyResult = 1;
                else
                    classifyResult = -1;
                end
                eachNetworkGuess(i,oo) = classifyResult;
            end
            
        end
        
        finalGuess = sum(eachNetworkGuess,2);
        for i = 1:length(finalGuess)
            if finalGuess(i)>0
                finalGuess(i) = 1;
            else
                finalGuess(i) = -1;
            end
        end
        
        correctVector = outcomesVal' - finalGuess;
        k = find(correctVector);
        percentCorrect = ((length(outcomesVal)-length(k))/length(outcomesVal)).*100;
        
        percentagesMatrix = [percentagesMatrix;percentCorrect];
        networksUsed{N,1} = networksChosen;
        
    end
end

[~,index2] = sort(percentagesMatrix(:,end),'descend');
percentagesMatrix = percentagesMatrix(index2,:);
networksUsed = networksUsed(index2);
%wrongVals = wrongVals(index2);


save(strcat('HyperParameterTests',filesep,'Full21Ensemble',featuresString,datestr(now),'.mat'),'percentagesMatrix','networksUsed')

profile viewer


