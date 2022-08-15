function [locs_start,locs_end]=amie_seg_test2(E,fs)
Th1=max(E)/4;
Th2=mean(E)/2;

d_Th=abs(Th1-Th2);
%% Adaptive threshold
Th_ad=0.5;
%Th_ad=mean(E)/3;
%d_Th=Th_ad/10;
Th_add = zeros(1,100); 
Th_add(1) = Th_ad;
cont = 2;
while Th_ad<Th2
    Th_ad = Th_ad + d_Th;
    Th_add(cont) = Th_ad;
    cont = cont + 1;
    if cont>100
        break
    end
end
Numb = length(find(Th_add>0));
Th_ad = Th_add(1:Numb);
%% Annotation (Silent Frame, Potential-breath Frame and Breath Frame)
Flag_S = zeros(Numb,length(E));
Flag_P = zeros(Numb,length(E));
Flag_B = zeros(Numb,length(E));
% Flagging frames
for i = 1:Numb
    %Silent frames
    Flag_S(i,:) = E<Th_ad(i);
    % Potential breath frames
    Flag_P(i,:) = E>Th_ad(i) & E<Th1;
    % Breath frames
    Flag_B(i,:) = E>Th1;
end



%POI=zeros(size(E)); 
% minimum length of segments is 31.25ms
thres=31.25/1000*fs; %currently around 10
delta_d=0.25*fs; %For paper experiments was 0.25 seconds
count=1;
for j=1:Numb 
    locs=find(Flag_S(j,:)); %silient frames
    for i=1:length(locs)-1
        if Flag_S(j,locs(i)+1)==1
            continue
        elseif Flag_B(j,locs(i)+1)==1
            continue
        elseif  (locs(i+1)-locs(i)-1)>=thres
            if sum(Flag_P(j,locs(i)+1:locs(i+1)-1))==(locs(i+1)-locs(i)-1)
                locs_start(count)=locs(i);
                locs_end(count)=locs(i+1); 
                count=count+1;
                %POI(locs(i))=1; %start
                %POI(locs(i+1))=-1; %end
            end
        else
            temp_seg=Flag_B(j,locs(i)+1:locs(i+1)-1);
            temp=find(temp_seg);
            if isempty(temp)
                continue
            elseif length(temp)==1
                locs_start(count)=locs(i);
                locs_end(count)=locs(i+1); 
                count=count+1;
                %POI(locs(i))=1; %start
                %POI(locs(i+1))=-1; %end
            elseif sum(diff(temp,1)==1)==(length(temp)-1)
                locs_start(count)=locs(i);
                locs_end(count)=locs(i+1); 
                count=count+1;
                %POI(locs(i))=1; %start
                %POI(locs(i+1))=-1; %end
            end
        end
    end
    
    %locs_start=find(POI==1);
    %locs_end=find(POI==-1);
    Dist_start=abs(diff(locs_start,1));
    Dist_end=abs(diff(locs_end,1));
    
    temp_start=(Dist_start<delta_d);
    temp_end=(Dist_end<delta_d);
    
    leave=sum(temp_start)+sum(temp_end);
    while leave>0
        loc1=find(temp_start,1); 
        loc2=find(temp_end,1); 
        
        if loc1<=loc2
            locs_start(loc1)=[];
            locs_end(loc1)=[];
        else
            locs_start(loc2+1)=[];
            locs_end(loc2+1)=[];
        end
        Dist_start=abs(diff(locs_start,1));
        Dist_end=abs(diff(locs_end,1));
        temp_start=(Dist_start<delta_d);
        temp_end=(Dist_end<delta_d);
        leave=sum(temp_start)+sum(temp_end);
    end
end
