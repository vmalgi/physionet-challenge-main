function y=hobalka(x,fs,fmin,fmax);
N=length(x);
    y=fft(x);
    if fmin>0
        imin=fix(fmin/(fs/N));
    else 
        imin=1;
        y(1)=y(1)/2;    %ss se neopakuje, nutno elimin. nas. 2 na konci
    end
    if fmax<fs/2
        imax=fix(fmax/(fs/N));
    else
        imax=fix(N/2);
    end
    
    hamwindow = hamming(imax-imin);
    hamsize = length(hamwindow);
    
    yy=zeros(size(y));
    
    istred=fix((imax+imin)/2);
    
    dolni=istred:imax;
    
    ld=length(dolni);
    
    yy(1:ld)=y(dolni).*hamwindow(floor(hamsize/2):hamsize);
    
    horni=imin:istred-1;
    lh=length(horni);
    
    yy(end-lh+1:end)=y(horni).*hamwindow(1:floor(hamsize/2));
    
    y=abs(ifft(yy))*2;
end

