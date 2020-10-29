%% population_UDS_analysis_closeloop.m
%% Collects summary data from UDS from anaesthetized recordings with closed-loop stimulation

%INPUT
%downsampled (200Hz) LFP recording in form 'xxxx_xH_rec1_downsampled.mat'
%US_indices.mat, generated my prelim_UDS_analysis.m

% important meta-data (recording hemisphere, baseline length) are hard-coded into
% script as 'allrecs' struct (see below) for each cohort

%OUTPUT
%(per cohort:)
%simple_incidence_all.mat %net upstate numbers
%incidence_all.mat %frequency of upstates;
%duration_all.mat %upstate duration
%auc_all.mat %area under the curve of upstates;
%gamma_all.mat %gamma power of upstates

%optional (per mouse): 'xxxx_recx_downstates.mat'
%downstate_struct.mat

%(c) Constance Holman
%%

%1. Cohort 1 =  Oxr1-Ai4D + CL
%2. Cohort 2 = Uchl1 + CL
%3. Cohort 3 = WT + CL


cohort = input('Enter which cohort you would like to analyze ');

switch cohort
    
    case 1 %Oxr1-Ai4D + CL
        basedir = 'your baseline directory';
        cd(basedir)
        allrecs = {007213 2 'LH' 10;
            007214 2 'RH' 5;...
            007216 2 'LH' 5;...
            007491 2 'LH' 10
            007489 2 'RH' 10};
        condition = cell2mat(allrecs(:,2)); %successful injections
        hemis = cell2mat(allrecs(:,3)); %checks which hemisphere recording occured in
        t_laser = cell2mat(allrecs(:,4)); %checks when stimulation started
        mousenrs = unique(cell2mat(allrecs(:,1)));
        fs = 200;
        
        n_recs =sum(cell2mat(allrecs(:,2)) == 2);
        
        
        for m = 1:length(allrecs)
            disp(allrecs(m,:))
            %go to correct directory
            mousestr = sprintf('DSC-00%d', allrecs{m,1});
            hemistr = hemis(m,:);
            dir = strcat('your baseline directory/Oxr1Ai4D_CL/',mousestr,'/',hemistr);
            cd(dir)
            datastr = strcat(mousestr,'_',hemistr(1:2),'_rec1_downsampled.mat');
            downstatestr = strcat(mousestr,hemistr(1:2),'_rec1_downstates.mat');
            load(datastr); load('US_indices.mat');   %   load(downstatestr);
            
            if m == 1
                baselinedownstates = downstate_struct.baseline;
                testdownstates = downstate_struct.test;
            else
                baselinedownstates = vertcat(baselinedownstates, downstate_struct.baseline);
                testdownstates = vertcat(testdownstates,downstate_struct.test);
            end
            
            
            n_recmins = ceil(length(dsdata)/fs/60); % find rec length in seconds
            
            mychan = 1; %unless latter decided otherwise...
            %TODO automate
            data = dsdata(mychan,:); %extracts one row
            
            baseline_min = t_laser(m);
            [incidence_temp, duration_temp, auc_temp, gamma_temp, downstate_struct, simple_incidence] = outsource_upstate_params(fs, baseline_min, indices, data);
            %
            incidence_all(m,:) = incidence_temp;
            duration_all(m,:) = duration_temp;
            auc_all(m,:) = auc_temp;
            gamma_all(m,:) = gamma_temp;
            simple_incidence_all(m,:) = simple_incidence;
            
            downstate_name = strcat(mousestr, hemistr, '_rec1_downstates.mat');
            save(downstate_name,'downstate_struct')
            
        end
        cd(basedir)
        % TODO
        
        
        simpleincidencename = 'Oxr1Ai4D_CL_simple_incidence';
        incidencename = 'Oxr1Ai4D_CL_incidence';
        durationname = 'Oxr1Ai4D_CL_duration';
        aucname = 'Oxr1Ai4D_CL_auc';
        gammaname = 'Oxr1Ai4D_CL_gamma';
        downstatename1 = 'Oxr1Ai4D_CL_downstates_baseline';
        downstatename2 = 'Oxr1Ai4D_CL_downstates_test';
%         save(downstatename1, 'baselinedownstates')
%         save(downstatename2, 'testdownstates')
%         %
%         save(simpleincidencename,'simple_incidence_all')
%         save(incidencename, 'incidence_all')
%         save(durationname, 'duration_all')
%         save(aucname, 'auc_all')
%         save(gammaname, 'gamma_all')
        %%
    case 2 %Uchl1 + CL
        
        basedir = 'your baseline directory';
        cd(basedir)
        
        allrecs = {010958 2 'RH' 10;010958 2 'LH' 10;...
            010959 2 'RH' 10;...
            010963 2 'RH' 10};
        
        condition = cell2mat(allrecs(:,2)); %good recs
        hemis = cell2mat(allrecs(:,3)); %checks which hemisphere recording occured in
        t_laser = cell2mat(allrecs(:,4)); %checks when stimulation started
        mousenrs = unique(cell2mat(allrecs(:,1)));
        fs = 200;
        
        n_recs =sum(cell2mat(allrecs(:,2)) == 2);
        
        
        
        for m = 1:length(allrecs)
            disp(allrecs(m,:))
            %go to correct directory
            mousestr = sprintf('DSC-0%d', allrecs{m,1});
            %hemistr = strcat(hemis(m,:), ' Spike Sorting');
            hemistr = hemis(m,:);
            dir = strcat('your baseline directory/Uchl1_CL/',mousestr,'/',hemistr);
            cd(dir)
            datastr = strcat(mousestr,'_',hemistr(1:2),'_rec1_downsampled.mat');
            downstatestr = strcat(mousestr,hemistr(1:2),'_rec1_downstates.mat');
            load(datastr); load('US_indices.mat'); %load(downstatestr);
            %
            %optional downstate analysis (files created in
            %'outsource_upstate_params'
%             if m == 1
%                 baselinedownstates = downstate_struct.baseline;
%                 testdownstates = downstate_struct.test;
%             else
%                 baselinedownstates = vertcat(baselinedownstates, downstate_struct.baseline);
%                 testdownstates = vertcat(testdownstates,downstate_struct.test);
%             end
            
            n_recmins = ceil(length(dsdata)/fs/60); % find rec length in seconds
            
            mychan = 1; %ch 9 is directly under optic fibre
            data = dsdata(mychan,:); %extracts one row
            baseline_min = t_laser(m);
            [incidence_temp, duration_temp, auc_temp, gamma_temp, downstate_struct, simple_incidence] = outsource_upstate_params(fs, baseline_min, indices, data);
            
            simple_incidence_all(m,:) = simple_incidence;
            incidence_all(m,:) = incidence_temp;
            duration_all(m,:) = duration_temp;
            auc_all(m,:) = auc_temp;
            gamma_all(m,:) = gamma_temp;
            downstate_name = strcat(mousestr, hemistr, '_rec1_downstates.mat');
            save(downstate_name,'downstate_struct')
        end
        cd('your baseline directory')

        simpleincidencename = 'Uchl1_CL_simple_incidence';
        incidencename = 'Uchl1_CL_incidence';
        durationname = 'Uchl1_CL_duration';
        aucname = 'Uchl1_CL_auc';
        gammaname = 'Uchl1_CL_gamma';
        downstatename1 = 'Uchl1_CL_downstates_baseline';
        downstatename2 = 'Uchl1_CL_downstates_test';
%         save(downstatename1, 'baselinedownstates')
%         save(downstatename2, 'testdownstates')
%         %
%         save(simpleincidencename, 'simple_incidence_all')
%         save(incidencename, 'incidence_all')
%         save(durationname, 'duration_all')
%         save(aucname, 'auc_all')
%         save(gammaname, 'gamma_all')
        
        downstate_name = strcat(mousestr, hemistr, '_rec1_downstates.mat');
        save(downstate_name,'downstate_struct')
        %%
        
    case 3 %WT cohort
        basedir = 'your baseline directory/Wildtype_CL';
        cd(basedir)
        allrecs = {44128 2 'RH' 10;...
            44129 2 'RH' 10; 044129 2 'LH' 10;...
            45790 2 'RH' 10;...
            45791 2 'RH' 10};
        condition = cell2mat(allrecs(:,2)); %good recs
        hemis = cell2mat(allrecs(:,3)); %checks which hemisphere recording occured in
        t_laser = cell2mat(allrecs(:,4)); %checks when stimulation started
        mousenrs = unique(cell2mat(allrecs(:,1)));
        fs = 200;
        
        n_recs =sum(cell2mat(allrecs(:,2)) == 2);
        
        indicence_all = NaN(n_recs,4);
        duration_all = NaN(n_recs,4);
        auc_all =NaN(n_recs,4);
        gamma_all = NaN(n_recs,4);
        
        
        for m = 1:size(allrecs,1)
            disp(allrecs(m,:))
            %go to correct directory
            mousestr = sprintf('SNA-0%d', allrecs{m,1});
            hemistr = hemis(m,:);
            dir = strcat('your baseline directory/Wildtype_CL/',mousestr,'/',hemistr);
            cd(dir)
            datastr = strcat(mousestr,'_',hemistr(1:2),'_rec1_downsampled.mat');
            %downstatestr = strcat(mousestr,hemistr(1:2),'_rec1_downstates.mat');

            load(datastr); load('US_indices.mat'); %load(downstatestr);
            
            %optional downstate analysis
%             if m == 1
%                 baselinedownstates = downstate_struct.baseline;
%                 testdownstates = downstate_struct.test;
%             else
%                 baselinedownstates = vertcat(baselinedownstates, downstate_struct.baseline);
%                 testdownstates = vertcat(testdownstates,downstate_struct.test);
%             end
            
            
            
            mychan = 1; %unless alter decided otherwise...
            data = dsdata(mychan,:); %extracts one row
            
            baseline_min = t_laser(m);
            
            [incidence_temp, duration_temp, auc_temp, gamma_temp, downstate_struct, simple_incidence] = outsource_upstate_params(fs, baseline_min, indices, data);
            
            simple_incidence_all(m,:) = simple_incidence;
            incidence_all(m,:) = incidence_temp;
            duration_all(m,:) = duration_temp;
            auc_all(m,:) = auc_temp;
            gamma_all(m,:) = gamma_temp;
            
%             downstate_name = strcat(mousestr, hemistr, '_rec1_downstates.mat');
%             save(downstate_name,'downstate_struct')
        end
        cd(basedir)
        %  % TODO
        simpleincidencename = 'WT_CL_simple_incidence';
        incidencename = 'WT_CL_incidence';
        durationname = 'WT_CL_duration';
        aucname = 'WT_CL_auc';
        gammaname = 'WT_CL_gamma';
        downstatename1 = 'WT_CL_downstates_baseline';
        downstatename2 = 'WT_CL_downstates_test';
%         save(downstatename1, 'baselinedownstates')
%         save(downstatename2, 'testdownstates')
%         %
%         
%         save(simpleincidencename, 'simple_incidence_all')
%         save(incidencename, 'incidence_all')
%         save(durationname, 'duration_all')
%         save(aucname, 'auc_all')
%         save(gammaname, 'gamma_all')
end

%% Plotting Summary Results

if cohort == 1
    groupstring = 'Oxr1-Ai4D CL';
elseif cohort == 2
    groupstring = 'Uchl1 CL'
elseif cohort == 3
    groupstring = 'WT_CL';
end

auc_all(isnan(auc_all)) = 0;
duration_all(isnan(duration_all)) = 0;
gamma_all(isnan(gamma_all)) = 0;
incidence_all(isnan(incidence_all)) = 0;

h1 = figure
subplot(1,4,1)
%
plot(incidence_all(:,[1,3])','k-o')
xtext={'Baseline' 'Stimulation'};
xlim([0 3])
set(gca,'xtick',[1:2],'xticklabel',xtext)
ylabel('Upstates per Minute (Hz)')
hold on
meanvals=mean(incidence_all);
plot(meanvals([1,3]),'--rs')
title('Incidence')
subplot(1,4,2)


%plot(incidence_all(:,[1,3])','k-o')
plot(simple_incidence_all','k-o')
xtext={'Baseline' 'Stimulation'};
xlim([0 3])
set(gca,'xtick',[1:2],'xticklabel',xtext)
ylabel('Average Upstates per Minute')
hold on
meanvals=mean(simple_incidence_all);
plot(meanvals,'--rs')
title('Incidence')
subplot(1,4,2)

yvec = auc_all(:,[1,3]);
xvec = repmat([1,2], 5,1);
err = auc_all(:,[2,4]);
%errorbar(xvec,yvec,err,'k-o')
plot(auc_all(:,[1,3])','k-o')
xtext={'Baseline' 'Stimulation'};
xlim([0 3])
set(gca,'xtick',[1:2],'xticklabel',xtext)
ylabel('uV2')
hold on
meanvals=mean(auc_all);
plot(meanvals([1,3]),'--rs')
title('Area Under The Curve')

subplot(1,4,3)

yvec = auc_all(:,[1,3]);
xvec = repmat([1,2], 5,1);
err = auc_all(:,[2,4]);
%errorbar(xvec,yvec,err,'k-o')
plot(gamma_all(:,[1,3])','k-o')
xtext={'Baseline' 'Stimulation'};
xlim([0 3])
set(gca,'xtick',[1:2],'xticklabel',xtext)
ylabel('Power (uV2)')
hold on
meanvals=mean(gamma_all);
plot(meanvals([1,3]),'--rs')
title('Mean Gamma Power in Upstates')

subplot(1,4,4)

plot(duration_all(:,[1,3])','k-o')
xtext={'Baseline' 'Stimulation'};
xlim([0 3])
set(gca,'xtick',[1:2],'xticklabel',xtext)
ylabel('Duration (ms)')
hold on
meanvals=mean(duration_all);
plot(meanvals([1,3]),'--rs')
title('Mean Duration of Upstates')

sgtitle(strcat('Summary: ',groupstring))

% h2 = figure
% cdfplot(baselinedownstates)
% hold on
% cdfplot(testdownstates)
% legend('Baseline','Stimulation')
% title(strcat(groupstring, ': Downstates'))
% xlabel('Duration (ms)')
% ylabel('Fraction of Time Lags')

% Summary Stats
summary_stats = struct('Incidence',signrank(incidence_all(:,1),incidence_all(:,3)));
summary_stats(1).AUC = signrank(auc_all(:,1),auc_all(:,3));
summary_stats(1).Duration = signrank(duration_all(:,1),duration_all(:,3));
summary_stats(1).Gamma = signrank(gamma_all(:,1),gamma_all(:,3));
% [h,p,ks2stats] =kstest2(baselinedownstates,testdownstates);
% summary_stats(1).DSlags = [h,p,ks2stats];

% summaryfigname = strcat(groupstring, ':_Summary');
% dsfigname = strcat(groupstring, ':_Downstates');
% statsname = strcat(groupstring, ': Statistics');

summaryfigname = 'your title here';
%TODO auomate with group names

statsname = strcat(groupstring, ': Statistics');
%%
savefig(h1,summaryfigname)
epscname = strcat(summaryfigname,'.epsc');
saveas(h1, epscname)
jpegname = strcat(summaryfigname,'.jpeg');
saveas(h1, jpegname)

%savefig(h2,dsfigname)

% saveas(h2, strcat(dsfigname,'.epsc'))
% saveas(h2, strcat(dsfigname,'.jpeg'))
save(statsname,'summary_stats')
