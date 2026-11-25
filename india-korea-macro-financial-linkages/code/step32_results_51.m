%% ================================================================
%  DESCRIPTIVE STATISTICS — for Results section 5.1
%  Reports mean, SD, min, max, skewness, kurtosis for each series,
%  in the form actually used in the models (differenced/growth rates).
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

% The variables as they enter the models (stationary forms)
vars = {'dln_nifty','dln_kospi','dln_india_ip','dln_korea_ip', ...
    'dln_usd_inr','dln_usd_kwr','dln_msci_world', ...
    'dln_india_cpi','dln_korea_cpi','dln_india_m3','dln_korea_m3', ...
    'd_indiacallrate','d_koreair','d_us3mtbills'};

fprintf('=== DESCRIPTIVE STATISTICS (variables in model form) ===\n\n');

results = {};
for i = 1:numel(vars)
    v = vars{i};
    x = data.(v);
    x = x(~isnan(x));
    n      = numel(x);
    mu     = mean(x);
    sd     = std(x);
    mn     = min(x);
    mx     = max(x);
    skw    = skewness(x);
    krt    = kurtosis(x);   % MATLAB reports raw kurtosis (normal = 3)
    results(end+1,:) = {v, n, mu, sd, mn, mx, skw, krt};
end

T = cell2table(results, 'VariableNames', ...
    {'Variable','N','Mean','SD','Min','Max','Skewness','Kurtosis'});
disp(T)

% Also report the two stock RETURN series and MSCI in annualised %
% terms for intuition (mean monthly log-diff * 12 * 100)
fprintf('\n--- For intuition: annualised mean return (%%) ---\n');
for v = {'dln_nifty','dln_kospi','dln_msci_world'}
    x = data.(v{1}); x = x(~isnan(x));
    fprintf('%-16s %.2f%% per year\n', v{1}, mean(x)*12*100);
end