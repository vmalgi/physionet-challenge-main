%% compute cwt sub-feature
function energytemp = cwtSubfeature(signal,timeN,freqN,Level,segind,ismean)
% input
% signal    input signal
% timeN     number of time interval 1,2,1,4
% freqN     number of frequency interval
% Level     number of resolution level
% segind    1 S1, 2 systolic 3 S2, 4 diastolic
% ismean    as claimed in cjy_cwt_feature

 freqintval = ComputeInt(Level,freqN);
 switch segind
     case 1 
         timeint = ComputeInt(size(signal,2),1*timeN);
     case 2
         timeint = ComputeInt(size(signal,2),2*timeN);
     case 3
         timeint = ComputeInt(size(signal,2),1*timeN);
     case 4
         timeint = ComputeInt(size(signal,2),4*timeN);
 end
 energytemp = zeros(freqN-1,length(timeint)-1);
 if(ismean)
      for i = 1:length(freqintval)-1
        for j = 1:length(timeint)-1
              energytemp(i,j) = sum(sum(abs(signal(freqintval(i):freqintval(i+1),timeint(j):timeint(j+1)))))/(timeint(j+1)-timeint(j));
        end
      end
 else
       for i = 1:length(freqintval)-1
        for j = 1:length(timeint)-1
              energytemp(i,j) = sum(sum(abs(signal(freqintval(i):freqintval(i+1),timeint(j):timeint(j+1)))));
        end
      end
 end 
end
  
    function intidx = ComputeInt(L,k) % separate L length into k intervals
              intidx = 1;
             for m = 1:k-1
                 intidx = [intidx round(L/k*m)];
             end
             intidx = [intidx L];
    end