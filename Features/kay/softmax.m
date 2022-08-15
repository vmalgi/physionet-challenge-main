function output = softmax(input)
    single = exp(input);
    output = single./(sum(single));
end