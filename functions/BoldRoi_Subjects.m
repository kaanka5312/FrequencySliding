function BOLD_matrix = BoldRoi_Subjects(subj)
%% =========================================================================
%  Script Name:    BoldRoi_Subjects.m
%  Author:         Kaan Keskin
%  Date Created:   2025-06-08
%  Last Modified:  2025-06-08
%
%  Description:
%  -------------------------------------------------------------------------
%
%  This script prepares the BOLD in the Glasser ROIs matrix for calculation
%  of further process such as Sliding Frequency, ACW and GSCORR
%
%  Inputs:
%    subj : subject folder name
%    
%
%  Outputs:
%    - 360 × T matrix per subject saved in BOLD_all_subjects cell array
%    - Optional 3D matrix: 360 × T × N subjects
%
%  Dependencies:
%    - MATLAB R2018 or later
%    - No external toolboxes required
%
%  -------------------------------------------------------------------------
%  Example Usage:
%    Run this script from the project root after setting up data folders.
%
%  -------------------------------------------------------------------------
%  Notes:
%    - Make sure .1D files contain numeric data only.
%    - Customize file and folder names if structure differs.
% =========================================================================
% Define ROI indices
roi_indices = [1:180, 1001:1180];

% Generate filenames
folder = [subj,'/Zscore_data_bandpass'];
filepaths = arrayfun(@(i) fullfile(folder, sprintf('BOLD_%d.1D', i)), roi_indices, 'UniformOutput', false);


% Load one file to get time length (assuming all same length)
BOLD_matrix = zeros(360, 158);  % 360 ROIs × T time points

% --- Load each file using filepaths ---
for i = 1:length(filepaths)
    if exist(filepaths{i}, 'file')
        data = load(filepaths{i}, '-ascii');
        BOLD_matrix(i, :) = data(:)';  % ensure row vector
    else
        warning('Missing file: %s', filepaths{i});
        BOLD_matrix(i, :) = nan(1, T);
    end
end

end