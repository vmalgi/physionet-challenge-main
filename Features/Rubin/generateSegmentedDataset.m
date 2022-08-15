% Generate segmented dataset
%
% Written by: Chengyu Liu, February 10 2016 chengyu.liu@emory.edu
%
% Last modified by: Jonathan Rubin, March 17 2016
%                   Jonathan.Rubin@parc.com
%

function [A, data] = generateSegmentedDataset(directory, fname)

	data_dir = [pwd filesep directory filesep];

	%% Add this directory to the MATLAB path.
	addpath(pwd)

    [A, data] = challenge([data_dir fname]);

