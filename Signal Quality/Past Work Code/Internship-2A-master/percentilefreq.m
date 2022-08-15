function p=percentilefreq(pxx,f,pi)
t=0;
for i=2:length(f)
t(i)=trapz(f(1:i),pxx(1:i));
end
t=100*t/max(t);
indx=find(t>=pi,1);
p=interp1(t(1:indx),f(1:indx),pi);

