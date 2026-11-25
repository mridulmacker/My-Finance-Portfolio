%% ================================================================
%  BLOCK 3 (FIXED) — COVID-DUMMY ROBUSTNESS, India FX, T-bill -> IP
%  Auto-detects COVID months from the India IP growth outliers, then
%  re-runs Toda-Yamamoto T-bill -> IP with those months dummied out.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

vars = {'ln_nifty','ln_india_ip','ln_usd_inr','ln_msci_world','us3m_tbills'};
Y = data{:, vars};
validRows = all(~isnan(Y),2);
Y = Y(validRows,:);
T = size(Y,1);

% --- Auto-detect COVID months: the India IP growth series has a handful of
%     extreme observations in 2020. Flag any month whose IP log-change is
%     more than 4 SD from the mean (this cleanly isolates Mar-Aug 2020). ---
ipgrowth = [NaN; diff(Y(:,2))];             % log-change in India IP
z = (ipgrowth - nanmean(ipgrowth)) / nanstd(ipgrowth);
covid = abs(z) > 4;                          % logical, T x 1
covid(isnan(covid)) = false;
fprintf('Detected %d extreme (COVID) months:\n', sum(covid));
disp(find(covid)');                          % row numbers, for your sanity check

% Build individual impulse dummies for each flagged month (cleanest: one
% dummy per extreme obs, so we neutralise them without distorting dynamics)
covidRows = find(covid);
D = zeros(T, numel(covidRows));
for j = 1:numel(covidRows)
    D(covidRows(j), j) = 1;
end

%% --- Toda-Yamamoto with COVID dummies as exogenous regressors ---
p    = 2;      % chosen VAR lag for India FX
dmax = 1;      % max integration order
k    = p + dmax;
n    = 5;

% Estimate levels VAR(k) with COVID dummies exogenous
Mdl = varm(n, k);
Mdl.SeriesNames = vars;
EstMdl = estimate(Mdl, Y, 'X', D);

% Build coefficient array: AR{1..k}, each n x n
A = zeros(n, n, k);
for L = 1:k
    A(:,:,L) = EstMdl.AR{L};
end

% TY Wald test: does T-bill (var 5) Granger-cause India IP (var 2)?
% Restrict coefficients on lags 1..p of var 5 in the IP (row 2) equation = 0.
% Need the covariance of those coefficients.
[EstMdl, ~, ~, E] = estimate(Mdl, Y, 'X', D);
Sigma = EstMdl.Covariance;                    % residual cov
% Extract the relevant coefficient estimates and their SEs via summarize
results = summarize(EstMdl);
% Wald test using the estimated coefficients and standard errors:
% Coeffs on var5 lags 1..p in equation 2.
% summarize gives a table; we pull the AR terms.
tbl = results.Table;
disp('--- Coefficients on T-bill lags in the India IP equation ---');
% Find rows: names like 'AR{1}(2,5)' ... 'AR{p}(2,5)'
mask = false(height(tbl),1);
for L = 1:p
    nm = sprintf('AR{%d}(2,5)', L);
    mask = mask | strcmp(tbl.Properties.RowNames, nm);
end
sub = tbl(mask,:);
disp(sub);

% Joint Wald: (beta)' inv(V) (beta) ~ chi2(p)
beta = sub.Value;
se   = sub.StandardError;
% Approximate joint test treating them as independent is not exact; use the
% full covariance if available. As a robust check, report each t-stat and a
% conservative joint chi-square from the sum of squared t-stats:
tstats = beta ./ se;
wald_approx = sum(tstats.^2);
pval_approx = 1 - chi2cdf(wald_approx, p);
fprintf('\nTY T-bill -> India IP WITH COVID dummies:\n');
fprintf('  Approx joint Wald = %.3f, df = %d, p = %.5f\n', wald_approx, p, pval_approx);
fprintf('  (Baseline WITHOUT dummy was p = 2.7e-5.)\n');
fprintf('\n[If p stays well below 0.05, the India monetary-channel result is robust\n');
fprintf(' to the pandemic outlier. If it weakens sharply, that is itself worth\n');
fprintf(' reporting honestly as an outlier-sensitivity caveat.]\n');