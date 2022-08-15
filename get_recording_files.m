function recording_files=get_recording_files(current_header)

recording_files={};

num_locations=strsplit(current_header{1},' ');
num_locations=str2double(num_locations{2});

for j=2:num_locations+1

    current_line=strsplit(current_header{j},' ');
    recording_files{j-1}=current_line{3};

end

end
