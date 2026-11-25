%% ================================================================
%  TIDY — Rebuild master_dataset_transformed.mat cleanly from scratch
%  ================================================================
clear; clc;
data = readtable('master_dataset.csv');
data.date = datetime(data.date);

log_vars = {'kospi','nifty','usd_inr','usd_kwr','msci_world', ...
    'india_ip','korea_ip','india_cpi','korea_cpi', ...
    'india_m3','korea_m3','india_wpi'};

for i = 1:numel(log_vars)
    v = log_vars{i};
    data.(['ln_' v]) = log(data.(v));
    data.(['dln_' v]) = [NaN; diff(data.(['ln_' v]))];
end

data.d_indiacallrate = [NaN; diff(data.india_call_rate)];
data.d_koreair        = [NaN; diff(data.korea_ir)];
data.d_us3mtbills     = [NaN; diff(data.us3m_tbills)];

fprintf('Clean rebuild: %d rows, %d columns\n', height(data), width(data));
disp(data.Properties.VariableNames')

save('master_dataset_transformed.mat', 'data');
fprintf('\nSaved clean master_dataset_transformed.mat\n');