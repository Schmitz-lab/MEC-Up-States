function [incidence_temp, duration_temp, auc_temp, gamma_temp, downstate_struct, simple_incidence...
    ] = outsource_upstate_params(fs, baseline_min, indices, data)

n_recmins = ceil(length(data)/fs/60); % find rec length in seconds

baseline_indices= [1,baseline_min*60*fs];
test_indices = [baseline_min*60*fs+1,baseline_min*60*fs+1 + baseline_min*60*fs+1]; %test period with laser on

baseline_upstates = find(indices(:,1) < baseline_indices(2)); % provides starting indices of all upstates in baseline
test_upstates = find(indices(:,1) > test_indices(1) & indices(:,1) < test_indices(2));

baseline_all = [indices(baseline_upstates,1),indices(baseline_upstates,2)];
test_all = [indices(test_upstates,1),indices(test_upstates,2)];

baseline_duration = baseline_all(:,2)-baseline_all(:,1);
test_duration = test_all(:,2)-test_all(:,1);

%%
auc_baseline = NaN(size(baseline_upstates));
auc_test = NaN(size(test_upstates));

if max(max(baseline_all)) > length(data)
    baseline_all = baseline_all([1:end-1],:)
end

for i = 1:length(auc_baseline)
    auc_baseline(i) = trapz(data(baseline_all(i,1):baseline_all(i,2))); %calculates simple area under the curve for upstates in baseline
end

for j = 1:length(auc_test)
    auc_test(j) = trapz(data(test_all(j,1):test_all(j,2)));
end

%% gamma etc
[AllMaxGamma] = gamma_UDS_Stockwell(data,fs);

gamma_baseline = NaN(size(baseline_upstates));
gamma_test = NaN(size(test_upstates));

for i = 1:length(gamma_baseline)
    gamma_baseline(i) = nanmean(AllMaxGamma(baseline_all(i,1):baseline_all(i,2))); % nanmean gamma power in upstates
end

for j = 1:length(gamma_test)
    gamma_test(j) = nanmean(AllMaxGamma(test_all(j,1):test_all(j,2)));
end

% calculate incidence


min_bins = 1:fs*60:baseline_min*fs*60;
incidence = nanmean((histcounts(baseline_all(:,1),min_bins))/60); %average incidence in hz
incidence_nanstd = nanstd((histcounts(baseline_all(:,1),min_bins))/60);
incidence_temp(1,1) = incidence;
incidence_temp(1,2) = incidence_nanstd;

simple_incidence(1,1) = length(baseline_all)/baseline_min;
simple_incidence(1,2) = length(test_all)/baseline_min;

test_min_bins = baseline_min*fs*60:fs*60:baseline_min*fs*60+baseline_min*fs*60;
incidence_test_mean = nanmean((histcounts(test_all(:,1),test_min_bins))/60); %average incidence in hz
incidence_test_nanstd = nanstd((histcounts(test_all(:,1),test_min_bins))/60);
incidence_temp(1,3) = incidence_test_mean;
incidence_temp(1,4) = incidence_test_nanstd;
%
%Optional plotting of upstates per min
% h=figure
% plot(histcounts(baseline_all(:,1),min_bins))
% hold on
% plot(histcounts(test_all(:,1),test_min_bins))
% legend('Baseline','Test')
% title('Upstate Incidence Per Min')
% savefig(h,'Upstate_Incidence')

duration = nanmean(baseline_duration);
duration_nanstd = nanstd(baseline_duration);
duration_temp(1,1) = duration;
duration_temp(1,2) = duration_nanstd;

duration_test_mean = nanmean(test_duration);
duration_test_nanstd =nanstd(test_duration);
duration_temp(1,3) = duration_test_mean;
duration_temp(1,4) = duration_test_nanstd;

auc = nanmean(auc_baseline);
auc_nanstd = nanstd(auc_baseline);
auc_temp(1,1) = auc;
auc_temp(1,2) = auc_nanstd;

auc_test_mean = nanmean(auc_test);
auc_test_nanstd =nanstd(auc_test);
auc_temp(1,3) = auc_test_mean;
auc_temp(1,4) = auc_test_nanstd;

gamma = nanmean(gamma_baseline);
gamma_nanstd = nanstd(gamma_baseline);
gamma_temp(1,1) = gamma;
gamma_temp(1,2) = gamma_nanstd;

gamma_test_mean = nanmean(gamma_test);
gamma_test_nanstd =nanstd(gamma_test);
gamma_temp(1,3) = gamma_test_mean;
gamma_temp(1,4) = gamma_test_nanstd;

% time between upstates
baseline_lags = baseline_all(:,2) - baseline_all(:,1);
test_lags = test_all(:,2) - test_all(:,1);
downstate_struct = struct('baseline',baseline_lags,'test',test_lags);

%optional plotting of lags between upstates
% figure
% subplot(1,2,1)
% plot(baseline_all(:,2) - baseline_all(:,1))
% hold on
% plot(test_all(:,2) - test_all(:,1))
% subplot(1,2,2)
% lags = [mean(baseline_all(:,2) - baseline_all(:,1)), mean(test_all(:,2) - test_all(:,1))];
% lags_std=[std(baseline_all(:,2) - baseline_all(:,1)), std(test_all(:,2) - test_all(:,1))];
% bar([1:2],lags)
% hold on
% errorbar(lags, lags_std)


