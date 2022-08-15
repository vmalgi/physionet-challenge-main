classdef constants
    %CONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Constant)
        minSysDuration = 175;
        maxSysDuration = 600; % ms
        
        minRRDuration = 300; %m�la by b�t krat��. ale nem��e b�t krat��. ne� min Systola
        maxRRDuration = 1800; 
        
        rrSysRatioLowLimit = 2.1; %limit pro korekci Sys/RR
        rrSysRatioHighLimit = 4.5; %limit pro korekci Sys/RR
        reliableRRsSTD = 3;
    end
    
    methods
    end
    
end

