
clearvars -except ROItimeseriesdata
idx = {'1' '2' '3' '4' '5' '6' '7' '8' '9' '10' '11' '12' '13' '14' '15' '16' '17' '18' '19' '20' '21' '22' '23' '24' '25' '26' '27' '28' '29' '30' '31' '32' '33' '34' '35' '36' '37' '38' '39' '40' '41' '42' '43' '44' '45' '46' '47' '48' '49' '50' '51' '52' '53' '54' '55' '56' '57' '58' '59' '60' '61' '62' '63' '64' '65' '66' '67' '68' '69' '70' '71' '72' '73' '74' '75' '76' '77' '78' '79' '80' '81' '82' '83' '84' '85' '86' '87' '88' '89' '90' '91' '92' '93'};
% ROI=8;
ROI=24;
% ROI=16;

all_pf = [];
pf_title = {};


for nROI = 1:ROI
valuelist = [];
title_list= {};
for sub = 1:numel(idx)
    %% get and set data
    %     index = ismember(rawdatasub(:,1),str2double(idx(sub))); %if your rawdata is double
    index = ismember(ROItimeseriesdata{:, 1},str2double(idx(sub))); %if your rawdata is table
    [row,col] = find(index);
    x = ROItimeseriesdata{row,2+nROI};
    varargin = 20;
    remove_time = 5;
    TR = 2;
    srate = 1 / TR;
    %% data setting
    [peak_freq,x_phase] = pf_cohen(x,srate,varargin,remove_time);
    

    % find the range of each phase
    idx_peak = find(x_phase >= -pi/12 & x_phase <= pi/12); %peak = −15~15 degree
    idx_peak = idx_peak - 1; % Adjust x_phase(2:end) to match peak_freq(1:end)
    idx_peak(idx_peak < 1) = []; % Ensure no index falls below 1

    idx_trough = find((x_phase >= 11*pi/12 & x_phase <= pi) | (x_phase >= -pi & x_phase <= -11*pi/12));%trough = 165~180 & -180 ~-165 degree
    idx_trough = idx_trough - 1;
    idx_trough(idx_trough < 1) = []; 

    idx_rise = find(x_phase >= -7*pi/12 & x_phase <= -5*pi/12);%rise =-105~-75 degree
    idx_rise = idx_rise - 1;
    idx_rise(idx_rise < 1) = [];

    idx_fall = find(x_phase >= 5*pi/12 & x_phase <= 7*pi/12);%fall = 75~105 degree
    idx_fall = idx_fall - 1;
    idx_fall(idx_fall < 1) = [];
  
    % find peak_freq coresponse to that phase and calculate mean
    PF_peak = mean(peak_freq(idx_peak));
    PF_trough = mean(peak_freq(idx_trough));
    PF_rise = mean(peak_freq(idx_rise));
    PF_fall = mean(peak_freq(idx_fall));
    PF_whole = mean(peak_freq);
    
    valuelist = [valuelist; PF_peak, PF_trough, PF_rise, PF_fall, PF_whole];

 
end
all_pf = [all_pf,valuelist];

title_list = [title_list, ...
    strcat(ROItimeseriesdata.Properties.VariableNames{2+nROI}, '_PF_peak'), ...
    strcat(ROItimeseriesdata.Properties.VariableNames{2+nROI}, '_PF_trough'), ...
    strcat(ROItimeseriesdata.Properties.VariableNames{2+nROI}, '_PF_rise'), ...
    strcat(ROItimeseriesdata.Properties.VariableNames{2+nROI}, '_PF_fall'), ...
    strcat(ROItimeseriesdata.Properties.VariableNames{2+nROI}, '_PF_whole')];


pf_title = [pf_title,title_list];

end

% Combine data and title as table
pf_data = array2table(all_pf, 'VariableNames',pf_title);

% Save the table to the Excel file
writetable(pf_data, 'D:\NTUH\drug\7.frequency_sliding\pf_data_slow2_20_5_new_PFwhole.xlsx');