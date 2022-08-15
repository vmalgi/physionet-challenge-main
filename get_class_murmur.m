function class=get_class_murmur(input_header)

current_header=get_header(input_header);

class=current_header(startsWith(current_header,'#Murmur'));
class=strsplit(class{1},':');
class=strtrim(class{2});

end