function sex=get_sex(header)

sex=header(startsWith(header,'#Sex'));
sex=strsplit(sex{1},':');
sex=strtrim(sex{2});

end