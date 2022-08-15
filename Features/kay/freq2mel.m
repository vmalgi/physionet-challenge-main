function mel = freq2mel(freq)

mel = 1125.*log(1 + freq./700);

end