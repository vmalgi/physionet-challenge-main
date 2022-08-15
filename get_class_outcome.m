function class=get_class_outcome(input_header)

current_header=get_header(input_header);

class=current_header(startsWith(current_header,'#Outcome'));
class=strsplit(class{1},':');
class=strtrim(class{2});

end