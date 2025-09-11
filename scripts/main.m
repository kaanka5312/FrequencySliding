%% Creating BOLD Matrix for Glasser ROI's for each subject 
%=+=+=+=+=+=+=+=+=+=+=+=+
%### B I P O L A R ######
%=+=+=+=+=+=+=+=+=+=+=+=+
% Linux 
%cd /media/kaanka5312/HD-B1/BPB_proc/

% MacOS
%cd /Volumes/HD-B1/BPB_proc/

% Windows 
cd E:/BPB_proc/

% List of subject directories (adjust pattern as needed)
subject_dirs = dir('sub-*.results');
subject_dirs_1 = {subject_dirs.name};
subject_dirs_full_1 = cellfun(@(s) fullfile('/Volumes/HD-B1/BPB_proc/', s), subject_dirs_1, 'UniformOutput', false);

% Ses-1
cd /Volumes/HD-B1/Thesis/SES-1_BIDS/derivatives/afni_proc/
%cd E:/Thesis/SES-1_BIDS/derivatives/afni_proc/ 
subject_dirs = dir('sub-*.results');
subject_dirs_2 = {subject_dirs.name};
subject_dirs_full_2 = cellfun(@(s) fullfile('/Volumes/HD-B1/Thesis/SES-1_BIDS/derivatives/afni_proc/', s), subject_dirs_2, 'UniformOutput', false);

% Merging 
bp_subject_dirs = [subject_dirs_1, subject_dirs_2];
bp_subject_dirs_full = [subject_dirs_full_1, subject_dirs_full_2];

% Apply BoldRoi_Subjects to each folder using cellfun
BOLD_all_subjects_BP = cellfun(@BoldRoi_Subjects, bp_subject_dirs_full, 'UniformOutput', false);

bp_subject_dirs = string(bp_subject_dirs);
writematrix(subject_dirs, 'bpb_strings.txt', 'Delimiter', 'tab');


save('/Users/kaankeskin/projects/FrequencySliding/data/output/BOLD_all_subjects_BP.mat',...
    "BOLD_all_subjects_BP",'-mat')

% Ses-2
cd E:/Thesis/SES-2_BIDS/derivatives/afni_proc/
subject_dirs = dir('sub-*.results');
subject_dirs_ses2 = {subject_dirs.name};
subject_dirs_full_ses2 = cellfun(@(s) fullfile('E:/Thesis/SES-2_BIDS/derivatives/afni_proc/', s), subject_dirs_ses2, 'UniformOutput', false);
% Apply BoldRoi_Subjects to each folder using cellfun
BOLD_ses2_BP = cellfun(@BoldRoi_Subjects, subject_dirs_full_ses2, 'UniformOutput', false);
save('C:/Users/kaank/OneDrive/Belgeler/GitHub/FrequencySliding/data/output/BOLD_ses2_BP.mat',...
    "BOLD_ses2_BP",'-mat')


subject_dirs = string(subject_dirs_ses2);
writematrix(subject_dirs, 'bpb_ses2_strings.txt', 'Delimiter', 'tab');

%=+=+=+=+=+=+=+=+=+=+=+=+
%### C O N T R O L ######
%=+=+=+=+=+=+=+=+=+=+=+=+
% Linux 
cd /media/kaanka5312/HD-B1/BIDS/derivatives/afni/

% Macos 
cd /Volumes/HD-B1/BIDS/derivatives/afni/
% List of subject directories (adjust pattern as needed)
% Concatenate the numeric ranges
subject_nums = [35:38, 40:53, 55:57, 59:70];

% Create folder names like 'sub-35.results', 'sub-36.results', ...
hc_subject_dirs = arrayfun(@(n) sprintf('sub-%d.results', n), subject_nums, 'UniformOutput', false);

% Apply BoldRoi_Subjects to each folder using cellfun
BOLD_all_subjects_HC = cellfun(@BoldRoi_Subjects, hc_subject_dirs, 'UniformOutput', false);
save('/home/kaanka5312/projects/FrequencySliding/data/output/BOLD_all_subjects_HC.mat',...
    "BOLD_all_subjects_HC",'-mat')

%% Bandpassing for slow4 and slow5 
%=+=+=+=+=+=+=+=+=+=+=+=+
%### B I P O L A R ######
%=+=+=+=+=+=+=+=+=+=+=+=+

% ##########################
% ### P E A K   F R E Q ####
% ##########################

load('./data/output/BOLD_all_subjects_BP.mat');

BOLD_filtered_slow4 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.027, 0.073, 1/3)', ...
    BOLD_all_subjects_BP, ...
    'UniformOutput', false);

BOLD_filtered_slow5 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.01, 0.027, 1/3)', ...
    BOLD_all_subjects_BP, ...
    'UniformOutput', false);

[peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{4}(:,4:153)',1/3);
phase_segmented_pf(x_phase, peak_freq)

PF_all_slow4 = cell(1, numel(BOLD_filtered_slow4));
for i=1:numel(BOLD_filtered_slow4)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow4{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_all_slow4{i} = PF_struct;
    PF_all_slow4{i}.subj_id = bp_subject_dirs{i};
end

PF_all_slow5 = cell(1, numel(BOLD_filtered_slow5));
for i=1:numel(BOLD_filtered_slow5)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_all_slow5{i} = PF_struct;
    PF_all_slow5{i}.subj_id = bp_subject_dirs{i};
end

save('./data/output/PF_BP.mat',"PF_all_slow4","PF_all_slow5",'-mat')

% Session 2
load('C:/Users/kaank/OneDrive/Belgeler/GitHub/FrequencySliding/data/output/BOLD_ses2_BP.mat')

BOLD_filtered_slow4 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.027, 0.073, 1/3)', ...
    BOLD_ses2_BP, ...
    'UniformOutput', false);

BOLD_filtered_slow5 = cellfun( ...
    @(x) bandpass_cheby1(x', 0.01, 0.027, 1/3)', ...
    BOLD_ses2_BP, ...
    'UniformOutput', false);

[peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{4}(:,4:153)',1/3);
phase_segmented_pf(x_phase, peak_freq)

PF_ses2_slow4 = cell(1, numel(BOLD_filtered_slow4));
for i=1:numel(BOLD_filtered_slow4)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow4{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_ses2_slow4{i} = PF_struct;
    PF_ses2_slow4{i}.subj_id = subject_dirs_ses2{i};
end

PF_ses2_slow5 = cell(1, numel(BOLD_filtered_slow5));
for i=1:numel(BOLD_filtered_slow5)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_ses2_slow5{i} = PF_struct;
    PF_ses2_slow5{i}.subj_id = subject_dirs_ses2{i};
end

save('./data/output/PF_ses2.mat',"PF_ses2_slow4","PF_ses2_slow5",'-mat')
% ##########################
% ### A C W  ####
% ##########################

%--- main script -------------------------------------------------------
fs = 1/3;                               % 3-s TR  → 0.333 Hz
BP_acw = cellfun(@(m) row_acw(m,fs), ...
                 BOLD_all_subjects_BP, ...
                 'UniformOutput',false);

% ##########################
% ### M F   ####
% ##########################
% ---------- main script ----------------------------------------------
fs        = 1/3;               % 3-s TR  → 0.333 Hz
freqrange = [0.01 0.1];        % typical fMRI band

BP_mf = cellfun(@(m) row_mf(m, fs, freqrange), ...
                BOLD_all_subjects_BP, ...
                'UniformOutput', false);


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
    PF_all_slow4{i}.subj_id = hc_subject_dirs{i};
end

PF_all_slow5 = cell(1, numel(BOLD_filtered_slow5));
for i=1:numel(BOLD_filtered_slow5)
    [peak_freq,x_phase] = pf_cohen(BOLD_filtered_slow5{i}(:,4:153)',1/3);
    PF_struct = phase_segmented_pf(x_phase, peak_freq);
    PF_all_slow5{i} = PF_struct;
    PF_all_slow5{i}.subj_id = hc_subject_dirs{i};
end

save('./data/output/PF_HC.mat',"PF_all_slow4","PF_all_slow5",'-mat')
