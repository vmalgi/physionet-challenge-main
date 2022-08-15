function current_recordings=load_recordings(input_directory,current_header)

recording_files=get_recording_files(current_header);

current_recordings={};

for j=1:length(recording_files)

    current_recordings{j}=audioread(fullfile(input_directory,strtrim(recording_files{j})));

end

end