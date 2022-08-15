
clear all;
close all;
clc

data_dir = [pwd filesep 'validation' filesep];

%% Add this directory to the MATLAB path.
addpath(pwd)

row = 0;


for p1=-0.5:0.1:0.5
    row = row +1;
    col = 0;
    
    for p2=p1:0.1:0.5
    col = col+1;
        
 fprintf('%nx%n','cell ',row,col);
    
fip=fopen('param1.m','wt');
 fprintf(fip,'%s\n','function res=param1');
 fprintf(fip,'%s%f;\n','res=',p1);
 fprintf(fip,'%s\n','end');
fclose(fip);

fip=fopen('param2.m','wt');
 fprintf(fip,'%s\n','function res=param2');
 fprintf(fip,'%s%f;\n','res=',p2);
 fprintf(fip,'%s\n','end');
fclose(fip);
        

test1 = param1();
test2 = param2();

%% Check for previous files before starting validation procedure
answers = dir(['answers.txt']);
if(~isempty(answers))
    while(1)
        display(['Found previous answer sheet file in: ' pwd])
%         cont = input('Delete it (Y/N/Q)?','s');
        cont = 'Y';
        if(strcmp(cont,'Y') || strcmp(cont,'N') || strcmp(cont,'Q'))
            if(strcmp(cont,'Q'))
                display('Exiting script!!')
                return;
            end
            break;
        end
    end
    if(strcmp(cont,'Y'))
        display('Removing previous answer sheet.')
        delete(answers.name);
    end
end

%% Load the list of records in the validation set.
fid = fopen([data_dir 'RECORDS'],'r');
if(fid ~= -1)
    RECLIST = textscan(fid,'%s');
else
    error(['Could not open ' data_dir 'RECORDS for scoring. Exiting...'])
end
fclose(fid);
RECORDS = RECLIST{1};

%% Running on the validation set and obtain the score results
classifyResult = zeros(length(RECORDS),1);
total_time     = 0;

fid=fopen('answers.txt','wt');
for i = 1:length(RECORDS)
    fname = RECORDS{i};
    tic;
    classifyResult(i) = challenge([data_dir fname]);

    % write the answer to answers.txt file
    fprintf(fid,'%s,%d\n',RECORDS{i},classifyResult(i));

    total_time = total_time+toc;
%     fprintf(['---Processed ' num2str(i) ' out of ' num2str(length(RECORDS)) ' records.\n'])
end
fclose(fid);


%% Scoring
modacc(row,col) = score2016ChallengeBrute;
save('bruteForce_limits_NO_GROUPS_small.mat','modacc');
fprintf('%s%n%n=%f1.4','cell ',row,col,modacc(row,col));
    end
    
end

fprintf('---DONE----');