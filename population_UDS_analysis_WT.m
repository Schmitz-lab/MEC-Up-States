%% population_UDS_analysis_WT

% used for examining properties of US in WT animals with no interventions

%INPUT
%downsampled (200Hz) LFP recording in form 'xxxx_xH_rec1_downsampled.mat'
%US_indices.mat, generated my prelim_UDS_analysis.m

% important meta-data (recording hemisphere, baseline length) are hard-coded into
% script as 'allrecs' struct (see below) for each cohort

%OUTPUT:
% matrix containing mean frequency and duration of upstates in 15 minute
% recordings

%(c) Constance Holman

allrecs = {58 1 'RH' 15 1 20 2;...
    62 1 'RH' 15 1 200 2;...
    64 1 'LH' 15 1 200 5;...
    65 1 'LH' 15 28 200 5;...
    70 1 'RH' 15 5 20 2;}

baseline_min = 15;

mousenrs = unique(cell2mat(allrecs(:,1)));
condition = cell2mat(allrecs(:,2)); %good recs
hemis = cell2mat(allrecs(:,3)); %checks which hemisphere recording occured in
baseline_end = cell2mat(allrecs(:,4)); %checks when baseline ended (all 15 in this case)
fs = 200;


% set up matrices to hold results
duration_all = NaN(size(allrecs,1), baseline_min);
frequency_all = NaN(size(allrecs,1), baseline_min);


% iterate through each mouse
for m = 1:size(allrecs,1)
    disp(allrecs(m,:))
    mousestr = sprintf('L%d', allrecs{m,1});
    hemis = cell2mat(allrecs(:,3)); %checks which hemisphere recording occured in
    hemistr = hemis(m,:); %modified 01.03.16 by CH to correspond to new directory system
    dir = ['your baseline directory/Wildtype/' mousestr '/' hemistr '/'];
    
    cd(dir)
    datastr = strcat(mousestr,'_',hemistr(1:2),'_rec1_downsampled_200Hz.mat');
    load(datastr); load('US_indices.mat');
    
    my_chan = dsdata(1,[1:baseline_min*60*fs]);
    last_index = length(indices(indices(:,1)<baseline_min*60*fs));
    
    baseline_us = indices([1:last_index],:);
    
    
    bin_edges = [1:fs*60:baseline_min*fs*60+1];
    incidence = histcounts(baseline_us(:,1),bin_edges)/60; %average incidence in hz
    frequency_all(m,:) = incidence;
    
    % sort US into minutes
    % find the duration of each US within a minute
    % find average of the durations and store them
    
    which_min = discretize(baseline_us(:,1),bin_edges);
    temp_durations = baseline_us(:,2) - baseline_us(:,1);
    mean_durations = NaN(1,baseline_min);
    
    for i =1:baseline_min
        mean_durations(i) = nanmean(temp_durations(which_min == i));
    end
    duration_all(m,:) = mean_durations;
end

cd('your baseline directory/Wildtype')

%TODO automatize names for saving variables?
name2save1 = 'WT_frequency_all.mat';
name2save2 = 'WT_duration_all.mat';
% save(name2save1, 'frequency_all')
% save(name2save2, 'duration_all')

