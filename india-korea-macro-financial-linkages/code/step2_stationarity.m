%% ================================================================
%  PART B — STATIONARITY TESTING (ADF) AND TRANSFORMATION
%  Replicates Code 1 (Stata) methodology:
%   - 12 lags throughout (Schwert 1989 rule of thumb for monthly data)
%   - Trend term for price/quantity levels; no trend for interest rates
%   - Differenced series tested with constant only (no trend), 12 lags
%  ================================================================

clear; clc;

%% Load data
data = readtable('master_dataset.csv');
data.date = datetime(data.date);

fprintf('Loaded %d observations, %d variables\n\n', height(data), width(data));

%% ---------- B3. Create log-differenced variables ----------
log_vars = {'kospi','nifty','usd_inr','usd_kwr','msci_world', ...
    'india_ip','korea_ip','india_cpi','korea_cpi', ...
    'india_m3','korea_m3','india_wpi'};

for i = 1:numel(log_vars)
    v = log_vars{i};
    ln_v = ['ln_' v];
    dln_v = ['dln_' v];
    data.(ln_v) = log(data.(v));
    data.(dln_v) = [NaN; diff(data.(ln_v))];
end

data.d_indiacallrate = [NaN; diff(data.india_call_rate)];
data.d_koreair        = [NaN; diff(data.korea_ir)];

%% ---------- B1. ADF tests on RAW LEVELS, 12 lags ----------
trended_levels = {'kospi','nifty','usd_inr','usd_kwr','msci_world', ...
    'india_ip','korea_ip','india_cpi','korea_cpi', ...
    'india_m3','korea_m3','india_wpi'};

nontrended_levels = {'india_call_rate','korea_ir'};

fprintf('=== B1. ADF on LEVELS (12 lags) ===\n\n');
results_levels = table();

for i = 1:numel(trended_levels)
    v = trended_levels{i};
    x = data.(v);
    x = x(~isnan(x));
    [h, pValue, stat] = adftest(x, 'lags', 12, 'model', 'TS');
    results_levels = [results_levels; {v, 12, "trend", stat, pValue, local_verdict(h)}];
end

for i = 1:numel(nontrended_levels)
    v = nontrended_levels{i};
    x = data.(v);
    x = x(~isnan(x));
    [h, pValue, stat] = adftest(x, 'lags', 12, 'model', 'ARD');
    results_levels = [results_levels; {v, 12, "none", stat, pValue, local_verdict(h)}];
end

results_levels.Properties.VariableNames = {'Variable','Lags','Trend','Zt_stat','pValue','Verdict'};
disp(results_levels)

%% ---------- B4. ADF tests on DIFFERENCED series, 12 lags ----------
diff_vars = {'dln_kospi','dln_nifty','dln_usd_inr','dln_usd_kwr','dln_msci_world', ...
    'dln_india_ip','dln_korea_ip','dln_india_cpi','dln_korea_cpi', ...
    'dln_india_m3','dln_korea_m3','dln_india_wpi', ...
    'd_indiacallrate','d_koreair'};

fprintf('\n=== B4. ADF on DIFFERENCED series (12 lags) ===\n\n');
results_diff = table();

for i = 1:numel(diff_vars)
    v = diff_vars{i};
    x = data.(v);
    x = x(~isnan(x));
    [h, pValue, stat] = adftest(x, 'lags', 12, 'model', 'ARD');
    results_diff = [results_diff; {v, 12, stat, pValue, local_verdict(h)}];
end

results_diff.Properties.VariableNames = {'Variable','Lags','Zt_stat','pValue','Verdict'};
disp(results_diff)

%% ---------- Cross-check: India M3 differenced at 4 lags ----------
fprintf('\n=== Cross-check: dln_india_m3 at 4 lags ===\n\n');
x = data.dln_india_m3;
x = x(~isnan(x));
[h, pValue, stat] = adftest(x, 'lags', 4, 'model', 'ARD');
fprintf('dln_india_m3, 4 lags: Z(t) = %.3f, p = %.3f -> %s\n', stat, pValue, local_verdict(h));

%% Save transformed dataset for later steps
save('master_dataset_transformed.mat', 'data');
fprintf('\nSaved master_dataset_transformed.mat\n');

%% ================================================================
function s = local_verdict(h)
if h == 1
    s = "Stationary";
else
    s = "Non-stationary";
end
end