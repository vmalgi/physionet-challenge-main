function features = get_features_2(header,recordings) %get_ECGLeads_features
%% Purpose
% Extract features from PCG recordings
%% Inputs
% header- header data
% recordings- The recordings
%% Output
% features- extracted features
%% Extract Header Information
age=get_age(header);
% neonate- birth to 27 days (0.5 months)
% infant- 28 days to 12 months (6 months)
% child- 2 to 11 years (6*12 months)
% adolescent- 12 years to 21 years (15*12 months)
% young adult- (20*12 months)
% nan- (0 months)

% sex=get_sex(header);
% % male (sex 0)
% % femae (sex 1)
% % nan (sex 2)
% 
% height=get_height(header);
% weight=get_weight(header);

recording_locations={'AV','MV','PV','TV','Phc'};
num_recording_locations=length(recording_locations);
%recording_features=zeros(num_recording_locations,4);

% Range of heart and respiratory heart for neonates

% neonate- birth to 27 days (0.5 months)
% infant- 28 days to 12 months (6 months)
% child- 2 to 11 years (6*12 months)
% adolescent- 12 years to 21 years (15*12 months)
% young adult- (20*12 months)
% nan- (0 months)
if strcmpi(age,'neonate')
    % birth to 27 days
    max_HR=200;
    min_HR=110;
    max_RR=80;
    min_RR=20;
elseif strcmpi(age,'infant')
    % 28 days to 12 months
    max_HR=200;
    min_HR=70;
    max_RR=60;
    min_RR=20;
elseif strcmpi(age,'child')
    % 2 years to 11 years
    max_HR=170;
    min_HR=60;
    max_RR=45;
    min_RR=15;
elseif strcmpi(age,'adolescent')
    % 12 years to 21 years
    max_HR=170;
    min_HR=40;
    max_RR=40;
    min_RR=10;
elseif strcmpi(age,'young adult')
    max_HR=130;
    min_HR=40;
    max_RR=40;
    min_RR=10;
else
    max_HR=160;
    min_HR=50;
    max_RR=80;
    min_RR=10;
end

locations=get_locations(header);
features_detailed=table; 
features=table; 
features_header=table;
SQI=table;
Fs=4000; 

options_hr.env='hilbert';
options_hr.autocorr='filtered';
options_hr.init_hr='autocorr_peak';
options_hr.systolic='yes';
options_hr.seg='springer';

load('Springer_B_matrix.mat', 'Springer_B_matrix');
load('Springer_pi_vector.mat', 'Springer_pi_vector');
load('Springer_total_obs_distribution.mat', 'Springer_total_obs_distribution');

options_hr.seg_fs=50;
options_hr.seg_pi_vector=Springer_pi_vector;
options_hr.seg_b_matrix=Springer_B_matrix;
options_hr.seg_total_obs_dist=Springer_total_obs_distribution;

first=0; 
for j=1:num_recording_locations
    fprintf('Recording: %d/%d \n',j,num_recording_locations)
    ind=strmatch(recording_locations{j},locations);
    if length(ind)==1
        audio_data=recordings{ind}; 
        %% Audio derived features
        % Mean
        % Example code
        features_header.age{j}=age;
%         features_header.sex{j}=sex;
%         features_header.height(j)=height; 
%         features_header.weight(j)=weight; 
        features_detailed.mean(j)=mean(abs(audio_data));
        
        temp = get_all_SQIs_2(audio_data,Fs,max_HR,min_HR,max_RR,min_RR);
        SQI(j,:)=splitvars(temp);
        
        %% Segmentation
        HR = get_hr_segmentation(audio_data,Fs,max_HR,min_HR,options_hr);
        %HR.seg_hr is the heart rate
        %HR.seg_states is the segmentation states of 1000Hz sampled version of
        %input heartsound 
        % There are 4 segmentation states as detailed below 
        % The S1 wave is identified by the integer 1.
        % The systolic period is identified by the integer 2.
        % The S2 wave is identified by the integer 3.
        % The diastolic period is identified by the integer 4.
        
        assigned_states=HR.seg_states{1}; 
        [PCG,Fs]=get_hr_preprocessing(audio_data,Fs,1000); 
        temp = get_heart_segmentation_features(assigned_states, PCG); 
        features(j,:)=splitvars(temp); 
        first=1; 
    else
        temp_features_detailed=table; 
        temp_SQI=table; 
        temp_features=table; 
        for i=1:length(ind)
            audio_data=recordings{ind(i)};
            temp_features_detailed.mean(i)=mean(abs(audio_data));
           
            temp = get_all_SQIs_2(audio_data,Fs,max_HR,min_HR,max_RR,min_RR);
            temp_SQI(i,:)= splitvars(temp);
            
            HR = get_hr_segmentation(audio_data,Fs,max_HR,min_HR,options_hr);
            assigned_states=HR.seg_states{1}; 
            [PCG,Fs]=get_hr_preprocessing(audio_data,Fs,1000); 
            temp = get_heart_segmentation_features(assigned_states, PCG); 
            temp_features(i,:)=splitvars(temp); 
        end
        if ~isempty(ind)
            features_detailed(j,:)=array2table(mean(temp_features_detailed{:,:}));
            SQI(j,:)=array2table(mean(temp_SQI{:,:}));
            features(j,:)=array2table(mean(temp_features{:,:}));
            features_header.age{j}=age;
%             features_header.sex{j}=sex;
%             features_header.height(j)=height;
%             features_header.weight(j)=weight;
            if first==0
                features_detailed.Properties.VariableNames=temp_features_detailed.Properties.VariableNames;
                SQI.Properties.VariableNames=temp_SQI.Properties.VariableNames;
                features.Properties.VariableNames=temp_features.Properties.VariableNames;
                first=1;
            end
        end
        if isempty(ind) && j==5
            features_detailed(6,:)=features_detailed(1,:);
            SQI(6,:)=SQI(1,:);
            features(6,:)=features(1,:);
            features_header(6,:)=features_header(1,:);
 
            features_detailed(6,:)=[];
            SQI(6,:)=[];
            features(6,:)=[];
            features_header(6,:)=[];
        end
        
    end
end


features=horzcat(features_header,features_detailed,SQI,features); 
end




































