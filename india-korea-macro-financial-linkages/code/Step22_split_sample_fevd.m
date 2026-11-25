%% ================================================================
%  SPLIT-SAMPLE FEVD CHECK — India C1 and Korea D1 primary specs
%  Compares FEVD across two clean, crisis-free sub-periods to test
%  whether the headline domestic:global finding is stable over time.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

function runSubsampleFEVD(data, vars, respVar, p, startDate, endDate, label)
mask = data.date >= startDate & data.date <= endDate;
Y = data{mask, vars};
validRows = all(~isnan(Y), 2);
Y = Y(validRows, :);
fprintf('\n--- %s (%s to %s), N=%d ---\n', label, datestr(startDate), datestr(endDate), size(Y,1));

if size(Y,1) < (p+1)*numel(vars)
    fprintf('Sample too small for reliable estimation at lag %d — skipping.\n', p);
    return;
end

Mdl = varm(numel(vars), p);
Mdl.SeriesNames = vars;
try
    [EstMdl, ~, ~, ~] = estimate(Mdl, Y);
catch ME
    fprintf('Estimation failed: %s\n', ME.message);
    return;
end

D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');
respIdx = find(strcmp(vars, respVar));
fevd_month12 = squeeze(D(12, :, respIdx));
fevdTable = table(vars', fevd_month12', 'VariableNames', {'ShockFrom','FEVD_Share'});
disp(fevdTable)
end

%% ================= INDIA C1 =================
fprintf('=====================================================\n');
fprintf('INDIA C1 (exchange rate spec) — split-sample FEVD check\n');
fprintf('Full-sample result (for reference): Nifty 11.8%%, MSCI World 10.7%%, own 85.6%%\n');
fprintf('=====================================================\n');

vars_ind = {'dln_nifty','dln_india_ip','dln_usd_inr','dln_msci_world'};
runSubsampleFEVD(data, vars_ind, 'dln_india_ip', 2, datetime(2000,1,1), datetime(2008,7,31), 'Pre-GFC');
runSubsampleFEVD(data, vars_ind, 'dln_india_ip', 2, datetime(2010,1,1), datetime(2019,12,31), 'Inter-crisis (post-GFC, pre-COVID)');

%% ================= KOREA D1 =================
fprintf('\n=====================================================\n');
fprintf('KOREA D1 (exchange rate spec) — split-sample FEVD check\n');
fprintf('Full-sample result (for reference): KOSPI 7.0%%, MSCI World 6.7%%, own 89.9%%\n');
fprintf('=====================================================\n');

vars_kor = {'dln_kospi','dln_korea_ip','dln_usd_kwr','dln_msci_world'};
runSubsampleFEVD(data, vars_kor, 'dln_korea_ip', 3, datetime(2000,1,1), datetime(2008,7,31), 'Pre-GFC');
runSubsampleFEVD(data, vars_kor, 'dln_korea_ip', 3, datetime(2010,1,1), datetime(2019,12,31), 'Inter-crisis (post-GFC, pre-COVID)');