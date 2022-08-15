function model = team_training_code(input_directory,output_directory) 
%% Purpose
% Train PCG classifiers and obtain the models
%% Inputs
% input_directory- folder containing the training set
% output_directory- folder to place the trained mde
%% Output
% model- the trained model

% Find text files
patient_files=dir(fullfile(input_directory,'*.txt'));
patient_files={patient_files.name};
patient_files=sort(patient_files); % To help debugging
num_patient_files=length(patient_files);

fprintf('Loading data for %d patients...\n', num_patient_files)

% Extract classes from data
classes_murmur={};
classes_outcome={};
for j=1:num_patient_files

    current_class_murmur=get_class_murmur(fullfile(input_directory,patient_files{j}));
    classes_murmur=unique([classes_murmur current_class_murmur]);

    current_class_outcome=get_class_outcome(fullfile(input_directory,patient_files{j}));
    classes_outcome=unique([classes_outcome current_class_outcome]);

end

classes_murmur=sort(classes_murmur);
num_classes_murmur=length(classes_murmur);

classes_outcome=sort(classes_outcome);
num_classes_outcome=length(classes_outcome);

% Extracting features and labels
features=table;
labels_murmur=categorical;
labels_outcome=categorical;

for j=1:num_patient_files
    fprintf('%d/%d \n',j,num_patient_files)

    current_header=get_header(fullfile(input_directory,patient_files{j}));
    current_recordings=load_recordings(input_directory,current_header);

    current_features=get_features_2(current_header, current_recordings);
    features((j-1)*5+1:j*5,:) = current_features;

    labels_murmur(j)=get_class_murmur(fullfile(input_directory,patient_files{j}));
    labels_outcome(j)=get_class_outcome(fullfile(input_directory,patient_files{j}));
end


%% Train model
all_labels_murmur=cell(5*length(labels_murmur),1);
all_labels_outcome=cell(5*length(labels_outcome),1);
all_files=zeros(5*length(labels_murmur),1);
for i=1:length(labels_murmur)
    all_labels_murmur{(i-1)*5+1}=char(labels_murmur(i));
    all_labels_murmur{(i-1)*5+2}=char(labels_murmur(i));
    all_labels_murmur{(i-1)*5+3}=char(labels_murmur(i));
    all_labels_murmur{(i-1)*5+4}=char(labels_murmur(i));
    all_labels_murmur{(i-1)*5+5}=char(labels_murmur(i));

    all_labels_outcome{(i-1)*5+1}=char(labels_outcome(i));
    all_labels_outcome{(i-1)*5+2}=char(labels_outcome(i));
    all_labels_outcome{(i-1)*5+3}=char(labels_outcome(i));
    all_labels_outcome{(i-1)*5+4}=char(labels_outcome(i));
    all_labels_outcome{(i-1)*5+5}=char(labels_outcome(i));
    
    all_files((i-1)*5+1)=str2double( erase(patient_files{i},'.txt'));
    all_files((i-1)*5+2)=str2double( erase(patient_files{i},'.txt'));
    all_files((i-1)*5+3)=str2double( erase(patient_files{i},'.txt'));
    all_files((i-1)*5+4)=str2double( erase(patient_files{i},'.txt'));
    all_files((i-1)*5+5)=str2double( erase(patient_files{i},'.txt'));
end
% remove empty rows
all_labels_murmur(features.mean==0)=[];
all_labels_outcome(features.mean==0)=[];
all_files(features.mean==0)=[];
features(features.mean==0,:)=[];

% normalise
[features(:,3:end),C,S] = normalize(features(:,3:end));


%% Murmur Model
[idx,~] = fscmrmr(features,all_labels_murmur,...
    'ClassNames',{'Absent','Present','Unknown'},...
    'Prior',[sum(strcmp(all_labels_murmur,'Absent')),...
        sum(strcmp(all_labels_murmur,'Present'))*5,...
        sum(strcmp(all_labels_murmur,'Unknown'))*3]);
% feature selection
top_features_murmur=idx(1:100);
top_train_features_murmur=features(:,top_features_murmur);


Cost_murmur=[0 1 1; 5 0 5; 3 3 0];
model_murmur = fitcensemble(top_train_features_murmur,all_labels_murmur,'Method','LPBoost',...
    'RatioToSmallest',[2.5 2.5 1],...
    'ClassNames',{'Absent','Present','Unknown'},'Cost',Cost_murmur,...
    'OptimizeHyperparameters',{'NumLearningCycles','LearnRate'},... 
    'HyperparameterOptimizationOptions',struct('Optimizer','bayesopt',...
    'Verbose',0,...
    'UseParallel',true,...
    'Repartition',false,...
    'ShowPlots',false));

%% Outcome Model
[idx,~] = fscchi2(features,all_labels_outcome,...
        'ClassNames',{'Abnormal','Normal'},...
        'Prior',[sum(strcmp(all_labels_outcome,'Abnormal'))*10,...
        sum(strcmp(all_labels_outcome,'Normal'))]);

% feature selection
top_features_outcome=idx(1:100);
top_train_features_outcome=features(:,top_features_outcome);

Cost_outcome=[0 10; 1 0];
model_outcome = fitcensemble(top_train_features_outcome,all_labels_outcome,'Method','LPBoost','ClassNames',{'Abnormal','Normal'},'Cost',Cost_outcome,...
    'OptimizeHyperparameters','none',... 
    'HyperparameterOptimizationOptions',struct('Optimizer','bayesopt',...
    'Verbose',0,...
    'UseParallel',true,...
    'Repartition',false,...
    'ShowPlots',false));

% save model
filename = fullfile(output_directory,'model.mat');
save(filename,'model_murmur','classes_murmur','model_outcome','classes_outcome',...
    'C','S','top_features_murmur','top_features_outcome','-v7.3');

end


