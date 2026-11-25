%% ================================================================
%  STEP 11 — INDIA C1 WITH US T-BILL ADDED (5-variable VAR)
%  First step: stationarity check on T-bill (never tested before).
%  Then: rebuild India C1 as Nifty, India IP, USD/INR, MSCI World, T-bill.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

%% ---- Stationarity check on US 3-month T-bill (previously untested) ----
fprintf('=== ADF: US 3-month T-bill, LEVEL (12 lags, no trend) ===\n');
x = data.us3m_tbills;
x = x(~isnan(x));
[h, pValue, stat] = adftest(x, 'lags', 12, 'model', 'ARD');
if h == 1
    verdict = "Stationary";
else
    verdict = "Non-stationary";
end
fprintf('Z(t) = %.4f, p = %.4f -> %s\n\n', stat, pValue, verdict);

% Create the differenced series
data.d_us3mtbills = [NaN; diff(data.us3m_tbills)];

fprintf('=== ADF: US 3-month T-bill, DIFFERENCED (12 lags) ===\n');
x = data.d_us3mtbills;
x = x(~isnan(x));
[h, pValue, stat] = adftest(x, 'lags', 12, 'model', 'ARD');
if h == 1
    verdict = "Stationary";
else
    verdict = "Non-stationary";
end
fprintf('Z(t) = %.4f, p = %.4f -> %s\n\n', stat, pValue, verdict);

save('master_dataset_transformed.mat', 'data');
fprintf('Updated master_dataset_transformed.mat with d_us3mtbills\n\n');

%% ================================================================
%  INDIA C1 + T-BILL (5 variables: Nifty, India IP, USD/INR, MSCI World, T-bill)
%  ================================================================
vars = {'dln_nifty','dln_india_ip','dln_usd_inr','dln_msci_world','d_us3mtbills'};
Y = data{:, vars};
validRows = all(~isnan(Y), 2);
Y = Y(validRows, :);
fprintf('Variables: %s\n', strjoin(vars, ', '));
fprintf('Full valid sample: %d observations\n\n', size(Y,1));

%% ---- Lag order selection ----
maxLag = 8;
nvar = numel(vars);
aicVals = nan(maxLag,1);
bicVals = nan(maxLag,1);
T_common = size(Y,1) - maxLag;

for p = 1:maxLag
    Mdl0 = varm(nvar, p);
    Ysub = Y(end-(T_common+p)+1:end, :);
    [EstMdl0, ~, logL] = estimate(Mdl0, Ysub);
    numParams = nvar*(nvar*p + 1);
    [aic, bic] = aicbic(logL, numParams, T_common);
    aicVals(p) = aic;
    bicVals(p) = bic;
end

lagTable = table((1:maxLag)', aicVals, bicVals, 'VariableNames', {'Lag','AIC','BIC'});
disp(lagTable)
[~, bestAIC] = min(aicVals);
fprintf('AIC selects lag = %d\n\n', bestAIC);

%% ---- Estimate VAR at AIC-selected lag ----
p = bestAIC;
Mdl = varm(nvar, p);
Mdl.SeriesNames = vars;
[EstMdl, EstSE, logL, Residuals] = estimate(Mdl, Y);
summarize(EstMdl)

%% ---- Stability check ----
n = nvar;
AR_all = [];
for i = 1:p
    AR_all = [AR_all, EstMdl.AR{i}];
end
if p > 1
    Companion = [AR_all; eye(n*(p-1)), zeros(n*(p-1), n)];
else
    Companion = AR_all;
end
E = eig(Companion);
fprintf('\nEigenvalue moduli (all should be < 1 for stability):\n');
disp(abs(E))
fprintf('Max modulus: %.4f\n', max(abs(E)));

%% ---- Residual autocorrelation ----
fprintf('\n=== Multivariate LM test for residual autocorrelation ===\n');
lmarTable = local_varlmar(Y, p, Residuals, 4);
disp(lmarTable)

%% ---- Granger causality ----
fprintf('\n=== Granger causality (leave-one-out, conditioned on full system) ===\n');
[h_gc, Summary_gc] = gctest(EstMdl, 'Display', false);
disp(Summary_gc)

%% ---- Generalized FEVD ----
fprintf('\n=== Generalized FEVD, response = India IP, horizon 12 ===\n');
D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');
respIdx = find(strcmp(vars, 'dln_india_ip'));
fevd_month12 = squeeze(D(12, :, respIdx));
fevdTable = table(vars', fevd_month12', 'VariableNames', {'ShockFrom','FEVD_Share'});
disp(fevdTable)

%% ---- Generalized IRF ----
irfResp = irf(EstMdl, 'NumObs', 12, 'Method', 'generalized');

save('india_c1_tbill_results.mat', 'EstMdl', 'Residuals', 'D', 'irfResp', 'Summary_gc', 'lmarTable');
fprintf('\nSaved india_c1_tbill_results.mat\n');


%% ================================================================
function tbl = local_varlmar(Y, p, U, mlagmax)
    [T, K] = size(Y);
    T_eff = T - p;
    X0 = ones(T_eff, 1);
    for j = 1:p
        X0 = [X0, Y(p+1-j:T-j, :)];
    end
    Sigma_hat = (U' * U) / T_eff;
    detSigma_hat = det(Sigma_hat);
    LM = nan(mlagmax,1); pval = nan(mlagmax,1); dfout = nan(mlagmax,1);
    for s = 1:mlagmax
        Us = zeros(T_eff, K);
        if s < T_eff
            Us(s+1:end, :) = U(1:end-s, :);
        end
        Z = [X0, Us];
        Beta = Z \ Y(p+1:T, :);
        E_s = Y(p+1:T, :) - Z * Beta;
        Sigma_tilde = (E_s' * E_s) / T_eff;
        d = K * size(Z, 2);
        LM(s) = (T_eff - d - 0.5) * log(detSigma_hat / det(Sigma_tilde));
        dfout(s) = K^2;
        pval(s) = 1 - chi2cdf(LM(s), K^2);
    end
    tbl = table((1:mlagmax)', LM, dfout, pval, 'VariableNames', {'Lag','LM_stat','df','pValue'});
end