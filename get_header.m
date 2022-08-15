function current_header=get_header(input_header)

current_header=fileread(input_header);
current_header=strsplit(current_header,'\n');

end