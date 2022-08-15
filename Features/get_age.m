function age_group=get_age(header)

age_group=header(startsWith(header,'#Age'));
age_group=strsplit(age_group{1},':');
age_group=strtrim(age_group{2});

end