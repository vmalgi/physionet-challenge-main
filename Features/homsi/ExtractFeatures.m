function ExtractFeatures()

    load('Springer_B_matrix.mat');
    load('Springer_pi_vector.mat');
    load('Springer_total_obs_distribution.mat');
    s = 'abcdef';
    for ns = 1:length(s) % Loop for each folder letter
       MyPath=['C:\Users\Masun\Documents\Physionet\training\training-',s(ns),'\'];
       FileName=[MyPath,'REFERENCE_withSQI.csv'];
%      FileName = 'C:\Users\Masun\Documents\Physionet\training\training-a\Reference.csv';
       fid = fopen(FileName);
       if ns<5
            out = textscan(fid,'%s%d%d%d%d','delimiter',',');
       else
            out = textscan(fid,'%s%d%d%d%d%d','delimiter',',');
       end;
       fclose(fid);
       springer_options   = default_Springer_HSMM_options;
       [L]=out{1};      % An array that has all files names
       [classes]=out{2};
       [Quality]=out{3};
       
       for i = 1:size(L) % Loop for each recording in a folder
           
            fn=strcat(MyPath,strtrim(L{i}));
            fn=strcat(fn,'.wav');
            MyClass=classes(i);
            MyQuality=Quality(i);
            MyClassQuality=strcat(num2str(MyClass),num2str(MyQuality));
            fprintf('%s',fn);

            [PCG, Fs1, nbits1] = wavread ([fn]);  % load data            
            PCG_resampled      = resample(PCG,springer_options.audio_Fs,Fs1);
            %% Running runSpringerSegmentationAlgorithm.m to obtain the assigned_states
            [assigned_states, heartRate] = runSpringerSegmentationAlgorithm(PCG_resampled, springer_options.audio_Fs, Springer_B_matrix, Springer_pi_vector, Springer_total_obs_distribution, false); % obtain the locations for S1, systole, s2 and diastole
    
            %% Running extractFeaturesFromHsIntervals.m to obtain the features for normal/abnormal heart sound classificaiton
            [features] = extractFeaturesFromHsIntervals(assigned_states,PCG_resampled, Fs1, heartRate);         
            allOneString = sprintf('%7.4f,' , features);
            allOneString = allOneString(1:end-1);
            allOneString1=strcat(strtrim(L{i}),',');
            allOneString1=strcat(allOneString1,allOneString);
            allOneString1=strcat(allOneString1,',');
            allOneString1=strcat(allOneString1,num2str(MyClass));
            allOneString1=strcat(allOneString1,',');
            allOneString1=strcat(allOneString1,num2str(MyQuality));
            allOneString1=strcat(allOneString1,',');
            allOneString1=strcat(allOneString1,num2str(MyClassQuality));
            
            if ns==1 && i==1
                %dlmwrite('myfile1.txt', pi, 'delimiter', '\t', 'precision', 16)
                dlmwrite('all-class-quality.csv',allOneString1, 'delimiter', '', 'precision', 16);          
            else
                dlmwrite('all-class-quality.csv',allOneString1,'-append', 'delimiter', '', 'precision', 16);               
            end;
            if ns==1 && i==1
                dlmwrite('a-class-quality.csv',allOneString1, 'delimiter', '', 'precision', 16);                   
            elseif ns==1 && i>1
                dlmwrite('a-class-quality.csv',allOneString1,'-append', 'delimiter', '', 'precision', 16);
            end;
            if ns==2 && i==1
                dlmwrite('b-class-quality.csv',allOneString1, 'delimiter', '', 'precision', 16);
            elseif ns==2 && i>1
                dlmwrite('b-class-quality.csv',allOneString1,'-append', 'delimiter', '', 'precision', 16);
            end;
            if ns==3 && i==1
                dlmwrite('c-class-quality.csv',allOneString1, 'delimiter', '', 'precision', 16);
            elseif ns==3 && i>1
                dlmwrite('c-class-quality.csv',allOneString1,'-append', 'delimiter', '', 'precision', 16);
            end;
            if ns==4 && i==1
                dlmwrite('d-class-quality.csv',allOneString1, 'delimiter', '', 'precision', 16);
            elseif ns==4 && i>1
                dlmwrite('d-class-quality.csv',allOneString1,'-append', 'delimiter', '', 'precision', 16);
            end;
            if ns==5 && i==1
                dlmwrite('e-class-quality.csv',allOneString1, 'delimiter', '', 'precision', 16);
            elseif ns==5 && i>1
                dlmwrite('e-class-quality.csv',allOneString1,'-append', 'delimiter', '', 'precision', 16);
            end;
            if ns==6 && i==1
                dlmwrite('f-class-quality.csv',allOneString1, 'delimiter', '', 'precision', 16);
            elseif ns==6 && i>1
                dlmwrite('f-class-quality.csv',allOneString1,'-append', 'delimiter', '', 'precision', 16);
            end;
            
                
            %%feat = getmswtfeat(PCG_Features,128,32,32)

    %         figure('Name', 'PCG features');

    %         t1 = (1:length(PCG_resampled))./springer_options.audio_Fs;
    %         plot(t1,PCG_resampled);
    %         hold on;
    %         t2 = (1:length(PCG_Features(:,4)))./featuresFs;
    %         plot(t2,PCG_Features(:,4));
            %%pause();
       end;
end

