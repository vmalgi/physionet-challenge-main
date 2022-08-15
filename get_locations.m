
function locations=get_locations(header)

num_locations=strsplit(header{1},' ');
num_locations=str2double(num_locations{2});

locations={};

for j=2:num_locations+1

    current_line=strsplit(header{j},' ');
    locations{j-1}=current_line{1};

end

end