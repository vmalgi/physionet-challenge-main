function class=get_class(input_header)
%% Purpose
% Get cardiac murmur class
%% Input
% input_header- location of header file
%% Output
% class- cardiac murmur present/absent/unknown
current_header=get_header(input_header);

class=current_header(startsWith(current_header,'#Murmur'));
class=strsplit(class{1},':');
class=strtrim(class{2});
end