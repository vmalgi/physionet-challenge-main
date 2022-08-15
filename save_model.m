function save_model(model_murmur,classes_murmur,model_outcome,classes_outcome,output_directory) %save_PCG_model
% Save results.
filename = fullfile(output_directory,'model.mat');
save(filename,'model_murmur','classes_murmur','model_outcome','classes_outcome','-v7.3');

disp('Done.')
end