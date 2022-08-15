for i=1:size(Lung_table,1)
     set=Lung_table.AnnotationSet(i);
     if set==1
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 1 Annotation/Lung Pool/';
     elseif set==2
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 2 Annotation/Lung Pool/';
     elseif set==3
         folder='/Users/ethangrooby/Downloads/Signal Quality/Annotation/Set 3 Annotation/Lung Pool/';
     end
     file=sprintf('%d.wav',Lung_table.index(i)); 
     [audio_data,Fs]=audioread([folder file]);
     
     RespiratoryRate=lung_sound_analysis_3(audio_data,Fs);
     duration=length(audio_data)/Fs;
     schmidt.pks(i,1)=length(RespiratoryRate.schmidt.overall_pks_locs);
     schmidt.hr(i,1)=RespiratoryRate.schmidt.overall_pks;
     schmidt.hr2(i,1)=RespiratoryRate.schmidt.overall_pks_2;
     Power.pks(i,1)=length(RespiratoryRate.Power.overall.peaks_locs);
     Power.hr(i,1)=RespiratoryRate.Power.overall.peaks;
     Power.hr2(i,1)=RespiratoryRate.Power.overall.peaks_2;
end


error=abs(schmidt.hr-BR_value);
 [error(Lung_qual>3),schmidt.hr(Lung_qual>3),BR_value(Lung_qual>3),find(Lung_qual>3)]
 median(error(Lung_qual==5))
 median(error(Lung_qual==4))
 
 error=abs(Power.hr-BR_value);
 [error(Lung_qual>3),Power.hr(Lung_qual>3),BR_value(Lung_qual>3),find(Lung_qual>3)]
 median(error(Lung_qual==5))
 median(error(Lung_qual==4))
 
 
  error=abs(Power.pks-BR_peaks);
 [error(Lung_qual>3),Power.pks(Lung_qual>3),BR_peaks(Lung_qual>3),find(Lung_qual>3)]
 median(error(Lung_qual==5))
 median(error(Lung_qual==4))
 
   error=abs(schmidt.pks-BR_peaks);
 [error(Lung_qual>3),schmidt.pks(Lung_qual>3),BR_peaks(Lung_qual>3),find(Lung_qual>3)]
 median(error(Lung_qual==5))
 median(error(Lung_qual==4))
 
 
 