function height=get_height(header)

height=header(startsWith(header,'#Height'));
height=strsplit(height{1},':');
height=strtrim(height{2});
height=str2double(height);

if isnan(height)
    height=0;
end

end