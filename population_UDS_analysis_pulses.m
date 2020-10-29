%% population_UDS_analysis_pulses

%INPUT
%downsampled (200Hz) LFP recording in form 'xxxx_xH_rec1_downsampled.mat'
%note chan#33 is from LED
%US_indices.mat, generated my prelim_UDS_analysis.m

% important meta-data (recording hemisphere, baseline length) are hard-coded into
% script as 'allrecs' struct (see below) for each cohort

%OUTPUT:
%'Oxr1Ai4D_pulses_simple_incidence.mat' column matrix of US incidence in
%LED off and LED on periods
%'Oxr1Ai4D_pulses_duration.mat' duration of upstates occurring in LED off
%and LED on periods

allrecs = {589 1 'RH' 15 32 20 2;... %nice unit
    585 1 'RH' 15 1 20 2;...
    585 1 'LH' 15 1 200 2;...
    584 1 'LH' 15 1 200 5;...
    720 1 'LH' 15 1 200 5;
    717 1 'LH' 15 1 200 5;
    717 1 'RH' 15 1 200 5};

condition = cell2mat(allrecs(:,2)); %successful injections
hemis = cell2mat(allrecs(:,3)); %checks which hemisphere recording occured in
t_led = cell2mat(allrecs(:,4)); %checks when stimulation started
mousenrs = unique(cell2mat(allrecs(:,1)));
fs = 200;


% iterate through mice
for m = 1:length(allrecs)
    disp(allrecs(m,:))
    %go to correct directory
    mousestr = sprintf('DSC-00E%d', allrecs{m,1});
    
    hemis = cell2mat(allrecs(:,3)); %checks which hemisphere recording occured in
    hemistr = hemis(m,:);
    dir = ['your baseline directory/Oxr1Ai4D_Pulses/' mousestr '/' hemistr];
    cd(dir)
    datastr = strcat(mousestr,'_',hemistr(1:2),'_rec1_downsampled.mat');
    
    load(datastr); load('US_indices.mat');
    
    %sanity plot for led pulses (optional)
    led_chan = dsdata(33,:);
    x = 1:length(led_chan);
    pulses = find(led_chan > 2.5);
    led_off= find(led_chan < 2.5);
    my_chan = dsdata(1,:);
    
    %%
    % cut out blocks of continuous led power (note/TODO: this could be done far
    % more simply with a threshold...)
    
    da=find(diff(pulses) ~= 1)+1; %find discontinuities in led "on" times
    da = [1, da, length(pulses)+1];
    
    for k = 1 : length(da)-1 %-2 compensates for also detecting "off" periods of led
        % First cell has the linear segments
        sa{k, 1} = pulses(da(k) : da(k+1)-1);
        % Second cell has the starting and stopping indices in the
        % original trace
        sa{k, 2} = [pulses(da(k)),pulses(da(k+1)-1)];
    end
    
    pulselength = 0.05 * fs; %50 ms pulse length in experiments
    n_pulses = length(pulses)/pulselength;
    
    %cuts random windows of 50ms when led is NOT on
    for p = 1:n_pulses
        index = randi(length(led_off)- pulselength);
        index2 = index + pulselength;
        off{p,1} = index:index2;
    end
    
    all_off = horzcat(off{1:end,1});
    
    
   %cycles through all detected upstates and tests whether they occur during a
    %pulse or not
    
   
    min_overlap = 0.05*fs; %minimum duration of upstates to be "counted" as occurring during pulse or not
    uds_count = 0;
    off_count = 0;
    duration_on=[];
    duration_off=[];
    
    for i = 1:length(indices)
        this_us = indices(i,1):indices(i,2);
        if sum(ismember(this_us,pulses))>= min_overlap %if more than 50ms of UDS are present in any given epoch??
            uds_count = uds_count + 1;
            duration_on(uds_count) = indices(i,2) - indices(i,1);
        end
        if sum(ismember(this_us,all_off))>= min_overlap %if more than 50ms of UDS are present in any given epoch??
            off_count = off_count + 1;
            duration_off  = indices(i,2) - indices(i,1);
        end
    end
    
    mean_duration_all(m,1) = mean(duration_off);
    mean_duration_all(m,2) = mean(duration_on);
    
    simple_incidence_all(m,1) = off_count/(length(all_off)/fs); %incidence in Hz
    simple_incidence_all(m,2) = uds_count/(length(pulses)/fs);
    
    
end % end of loop going through all mice
cd('your baseline directory/Oxr1Ai4D_Pulses')


simpleincidencename = 'Oxr1Ai4D_pulses_simple_incidence';
durationname = 'Oxr1Ai4D_pulses_duration';

%
save(simpleincidencename,'simple_incidence_all')
save(durationname, 'mean_duration_all')


