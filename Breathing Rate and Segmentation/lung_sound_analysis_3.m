function RespiratoryRate=lung_sound_analysis_3(x,fs)

x(x==0)=min(abs(x(x~=0)));  %dealing with padded zeros 
%x(x==0)=eps; %dealing with padded zeros 
x = resample(x,2000,fs);
fs=2000; 

[Power,~]=power_spectrum(x,fs);

duration=length(x)/fs;

RespiratoryRate.schmidt=computer_RR_schmidt(x,fs);
RespiratoryRate.Power=computer_RR(Power{2},10,duration);
end

function RespiratoryRate=computer_RR_schmidt(audio_data,fs)
[RespiratoryRate.overall_pks,heartRate_pks_locs] = getHeartRateSchmidt(audio_data, fs);
 
RespiratoryRate.overall_pks=RespiratoryRate.overall_pks;
RespiratoryRate.overall_pks_locs=heartRate_pks_locs;
RespiratoryRate.overall_pks_2=60*2000/mean(diff(heartRate_pks_locs));
end
function RespiratoryRate=computer_RR(audio_data,fs,duration)
%audio_data=normalize(audio_data);
audio_data=audio_data/max(audio_data); 
[~, locs] = findpeaks(audio_data,'MinPeakDistance' ,fs/(120/60),'MinPeakHeight',max(audio_data)/30,'MinPeakProminence', 0.03);

RespiratoryRate.overall.peaks_locs=locs;
RespiratoryRate.overall.peaks=fs*length(locs)/length(audio_data)*60;
RespiratoryRate.overall.peaks_2=60*fs/mean(diff(locs));
end

function [overall,t]=power_spectrum(audio_data,Fs)
%Computerised acoustical respiratory phase detection without airflow
%measurement paper Z. K. Moussavi et al.
%Bandpass 50-2500Hz
audio_data = butterworth_low_pass_filter(audio_data,2,min(2500,Fs/2-1),Fs);
audio_data = butterworth_high_pass_filter(audio_data,3,50,Fs);
% Hanning window, 50% overlap, 200ms segments
%200ms segment chosen as approximate duration of one breath cycle
window_length=200/1000*Fs; %200 ms window
overlap_length=window_length*0.5; %50% overlap
[s,f,t]=stft(audio_data,Fs,'Window',hann(window_length,'periodic'),'OverlapLength',overlap_length,'Centered',1==0);
% 150-450Hz gives the biggest difference in inspiration and expiration
% Another paper Breath Analysis of Respiratory Flow using Tracheal Sounds
% Saiful Huq et al. stated 300-450Hz had the biggest difference.
% Detected inspiratory peaks for analysis.
% They also considered 70-300, 600-800, 800-1000 and 1000-1200
low=[150 300 450 600];
high=[300 450 600 1200];
overall=cell(4,1);
for i=1:4
    begin=find(f>=low(i),1,'first');
    fin=find(f<=high(i),1,'last');
    temp=s(begin:fin,:);
    overall{i}=mean(abs(temp).^2);
end
end

function [heartRate_pks,heartRate_pks_locs] = getHeartRateSchmidt(audio_data, Fs)
%% 25-400Hz 4th order Butterworth band pass
%changed to 150-1000
%audio_data = butterworth_low_pass_filter(audio_data,2,1000,Fs, false);
audio_data = butterworth_high_pass_filter(audio_data,2,150,Fs);

%% Spike removal from the original paper:
%audio_data = schmidt_spike_removal(audio_data,Fs); %seems to have no effect

%% Find the homomorphic envelope
homomorphic_envelope = Homomorphic_Envelope_with_Hilbert(audio_data, Fs);

%% Find the autocorrelation:
%y=homomorphic_envelope-mean(homomorphic_envelope);
y=normalize(homomorphic_envelope);
[~, locs] = findpeaks(y,'MinPeakDistance' ,Fs/(120/60),'MinPeakHeight',max(y)/30,'MinPeakProminence', 0.03);
duration=length(y)/Fs;
heartRate_pks=length(locs)/duration*60;
heartRate_pks_locs=locs;
end
