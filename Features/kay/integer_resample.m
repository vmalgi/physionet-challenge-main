function [resampled_assigned_states] = integer_resample(assigned_states,downsampleFs,Fs)

resampled_assigned_states = resample(assigned_states,downsampleFs,Fs);
%resample does not quite give integers
resampled_assigned_states = int64(resampled_assigned_states);
%remove zeros with 1
resampled_assigned_states( resampled_assigned_states==0 )=1;
resampled_assigned_states( resampled_assigned_states==5 )=4;
%Problems when going from 4 to 1 when resampling
index = find((diff(resampled_assigned_states)==-2));
resampled_assigned_states( index+1 ) = 1;
index = find((diff(resampled_assigned_states)==-1));
resampled_assigned_states( index+1) = 1;

end