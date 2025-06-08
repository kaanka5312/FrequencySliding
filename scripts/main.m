%% Creating BOLD Matrix for Glasser ROI's for each subject 
%=+=+=+=+=+=+=+=+=+=+=+=+
%### B I P O L A R ######
%=+=+=+=+=+=+=+=+=+=+=+=+
% Linux 
cd /media/kaanka5312/HD-B1/BPB_proc/
% List of subject directories (adjust pattern as needed)
subject_dirs = dir('sub-*.results');
subject_dirs = {subject_dirs.name};

% Apply BoldRoi_Subjects to each folder using cellfun
BOLD_all_subjects_BP = cellfun(@BoldRoi_Subjects, subject_dirs, 'UniformOutput', false);
save('/home/kaanka5312/projects/FrequencySliding/data/output/BOLD_all_subjects_BP.mat',...
    "BOLD_all_subjects_BP",'-mat')

%=+=+=+=+=+=+=+=+=+=+=+=+
%### C O N T R O L ######
%=+=+=+=+=+=+=+=+=+=+=+=+
cd /media/kaanka5312/HD-B1/BIDS/derivatives/afni/
% List of subject directories (adjust pattern as needed)
% Concatenate the numeric ranges
subject_nums = [35:38, 40:53, 55:57, 59:70];

% Create folder names like 'sub-35.results', 'sub-36.results', ...
subject_dirs = arrayfun(@(n) sprintf('sub-%d.results', n), subject_nums, 'UniformOutput', false);

% Apply BoldRoi_Subjects to each folder using cellfun
BOLD_all_subjects_HC = cellfun(@BoldRoi_Subjects, subject_dirs, 'UniformOutput', false);
save('/home/kaanka5312/projects/FrequencySliding/data/output/BOLD_all_subjects_HC.mat',...
    "BOLD_all_subjects_HC",'-mat')

%% Bandpassing for slow4 and slow5 
%=+=+=+=+=+=+=+=+=+=+=+=+
%### B I P O L A R ######
%=+=+=+=+=+=+=+=+=+=+=+=+
BOLD_all_subjects=load('./data/output/BOLD_all_subjects_BP.mat');

BOLD_filtered_slow4 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.027, 0.073, 1/3)', ...
    BOLD_all_subjects, ...
    'UniformOutput', false);

BOLD_filtered_slow5 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.01, 0.027, 1/3)', ...
    BOLD_all_subjects, ...
    'UniformOutput', false);

[peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{4}(:,4:153)',1/3);
phase_segmented_pf(x_phase, peak_freq)

PF_all_slow4 = cell(1, numel(BOLD_filtered_slow4));
for i=1:numel(BOLD_filtered_slow4)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow4{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_all_slow4{i} = PF_struct;
end

PF_all_slow5 = cell(1, numel(BOLD_filtered_slow5));
for i=1:numel(BOLD_filtered_slow5)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_all_slow5{i} = PF_struct;
end

save('./data/output/PF_BP.mat',"PF_all_slow4","PF_all_slow5",'-mat')

%=+=+=+=+=+=+=+=+=+=+=+=+
%### C O N T R O L ######
%=+=+=+=+=+=+=+=+=+=+=+=+

load('./data/output/BOLD_all_subjects_HC.mat');
BOLD_all_subjects = BOLD_all_subjects_HC;

BOLD_filtered_slow4 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.027, 0.073, 1/3)', ...
    BOLD_all_subjects, ...
    'UniformOutput', false);

BOLD_filtered_slow5 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.01, 0.027, 1/3)', ...
    BOLD_all_subjects, ...
    'UniformOutput', false);

[peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{4}(:,4:153)',1/3);
phase_segmented_pf(x_phase, peak_freq)

PF_all_slow4 = cell(1, numel(BOLD_filtered_slow4));
for i=1:numel(BOLD_filtered_slow4)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow4{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_all_slow4{i} = PF_struct;
end

PF_all_slow5 = cell(1, numel(BOLD_filtered_slow5));
for i=1:numel(BOLD_filtered_slow5)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_all_slow5{i} = PF_struct;
end

save('./data/output/PF_HC.mat',"PF_all_slow4","PF_all_slow5",'-mat')
