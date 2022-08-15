%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Purpose: Run trained PCG classifier and obtain classifier outputs
% Inputs:
% 1. header
% 2. recordings
% 3. trained model
%
% Outputs:
% 1. score
% 2. label
% 3. classes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [score, label, classes] = team_testing_code(header, recordings, loaded_model)
classes_murmur=loaded_model.classes_murmur;
classes_outcome=loaded_model.classes_outcome;

features=get_features_2(header,recordings);
features(features.mean==0,:)=[];

% normalise
features(:,3:end) = normalize(features(:,3:end),'center',loaded_model.C,'scale',loaded_model.S);

%% Murmur
% feature selection
top_test_features_murmur=features(:,loaded_model.top_features_murmur);

[r,~]=size(top_test_features_murmur);
for i=1:r
    [label_tmp(i),score_murmur(i,:)]=loaded_model.model_murmur.predict(top_test_features_murmur(i,:));
end

guess=[sum(strcmp(label_tmp,'Present')) sum(strcmp(label_tmp,'Unknown')) sum(strcmp(label_tmp,'Absent'))];

final_guess=find(guess==max(guess),1);

label_murmur=zeros(1,length(classes_murmur));
if final_guess==1
    label_murmur(2)=1;
elseif final_guess==2
    label_murmur(3)=1;
else
    label_murmur(1)=1;
end

mean_score_murmur=mean(score_murmur,1); 


%% Outcome
% feature selection
top_test_features_outcome=features(:,loaded_model.top_features_outcome);

[r,~]=size(top_test_features_outcome);
for i=1:r
    [label_tmp(i),score_outcome(i,:)]=loaded_model.model_outcome.predict(top_test_features_outcome(i,:));
end

guess=[sum(strcmp(label_tmp,'Abnormal')) sum(strcmp(label_tmp,'Normal'))];


if guess(1)>0
    final_guess=1;
else
    final_guess=2;
end


label_outcome=zeros(1,length(classes_outcome));
if final_guess==1
    label_outcome(1)=1;
else
    label_outcome(2)=1;
end

mean_score_outcome=mean(score_outcome,1); 

classes=[classes_murmur classes_outcome];
score=[mean_score_murmur mean_score_outcome];
label=[label_murmur label_outcome];

end
