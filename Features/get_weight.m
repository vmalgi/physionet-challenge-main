function weight=get_weight(header)

weight=header(startsWith(header,'#Weight'));
weight=strsplit(weight{1},':');
weight=strtrim(weight{2});
weight=str2double(weight);

if isnan(weight)
    weight=0;
end

end