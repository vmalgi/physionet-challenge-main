function freq = mel2freq(mel)
    freq = 700.*(exp(mel./1125) - 1);
end