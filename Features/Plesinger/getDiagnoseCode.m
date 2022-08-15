function [ code ] = getDiagnoseCode(diagname)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

diagname = strtrim(diagname);
diagname = upper(diagname);

code = 0;

 if strcmp(diagname,'NORMAL') | strcmp(diagname,'CONTROL') | strcmp(diagname,'CONTROLS')
     code = -1;
 end
 
  if strcmp(diagname,'AD') 
     code = 1;
  end
 
   if strcmp(diagname,'AS') 
     code = 2;
   end
 
    if strcmp(diagname,'BENIGN') 
     code = 3;
    end
 
     if strcmp(diagname,'CAD') 
     code = 4;
     end
 
      if strcmp(diagname,'MPC') 
     code = 5;
 end
 
  if strcmp(diagname,'MR') 
     code = 6;
  end
 
   if strcmp(diagname,'MVP') 
     code = 7;
   end
 
    if strcmp(diagname,'PATHOLOGIC') 
     code = 8;
    end
 
  

end

