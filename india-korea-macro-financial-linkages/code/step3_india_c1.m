%% ================================================================
%  PART C1 — INDIA PRIMARY VAR (exchange rate spec)
%  Uses GENERALIZED impulse responses/FEVD (Pesaran-Shin, 1998).
%  Includes the CORRECT multivariate LM test for residual
%  autocorrelation (Johansen 1995), matching Stata's varlmar exactly:
%  LM_s = (T - d - 0.5) * ln(|Sigma_hat| / |Sigma_tilde_s|), df = K^2
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

vars = {'dln_nifty','dln_india_ip','dln_usd_inr','dln_msci_world'};
Y = data{:, vars};
validRows = all(~isnan(Y), 2);
Y = Y(validRows, :);
fprintf('Full valid sample: %d observations\n\n', size(Y,1));

%% ---- Lag order selection (varsoc equivalent) ----
maxLag = 8;
aicVals = nan(maxLag,1);
bicVals = nan(maxLag,1);
T_common = size(Y,1) - maxLag;

for p = 1:maxLag
    Mdl0 = varm(4, p);
    Ysub = Y(end-(T_common+p)+1:end, :);
    [EstMdl0, ~, logL] = estimate(Mdl0, Ysub);
    numParams = 4*(4*p + 1);
    [aic, bic] = aicbic(logL, numParams, T_common);
    aicVals(p) = aic;
    bicVals(p) = bic;
end

lagTable = table((1:maxLag)', aicVals, bicVals, 'VariableNames', {'Lag','AIC','BIC'});
disp(lagTable)
[~, bestAIC] = min(aicVals);
fprintf('AIC selects lag = %d\n\n', bestAIC);

%% ---- Estimate final VAR (lag 2, per Stata spec) ----
p = 2;
Mdl = varm(4, p);
Mdl.SeriesNames = vars;
[EstMdl, EstSE, logL, Residuals] = estimate(Mdl, Y);
summarize(EstMdl)

%% ---- Stability check (companion matrix) ----
n = 4;
AR1 = EstMdl.AR{1};
AR2 = EstMdl.AR{2};
Companion = [AR1, AR2; eye(n), zeros(n,n)];
E = eig(Companion);
fprintf('\nEigenvalue moduli (all should be < 1 for stability):\n');
disp(abs(E))
fprintf('Max modulus: %.4f\n', max(abs(E)));

%% ---- Residual autocorrelation: PROPER multivariate LM test (varlmar) ----
fprintf('\n=== Multivariate LM test for residual autocorrelation (Johansen 1995) ===\n');
fprintf('(exact replication of Stata''s varlmar formula, testing lags 1-4)\n\n');
lmarTable = local_varlmar(Y, p, Residuals, 4);
disp(lmarTable)

%% ---- Granger causality: SYSTEM-CONDITIONED (matches Stata's vargranger) ----
fprintf('\n=== Granger causality (leave-one-out, conditioned on full system) ===\n');
[h_gc, Summary_gc] = gctest(EstMdl, 'Display', false);
disp(Summary_gc)

%% ---- GENERALIZED FEVD (Pesaran-Shin, ordering-invariant) ----
fprintf('\n=== Generalized FEVD, response = India IP, horizon 12 ===\n');
D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');
respIdx = find(strcmp(vars, 'dln_india_ip'));
fevd_month12 = squeeze(D(12, :, respIdx));

fevdTable = table(vars', fevd_month12', 'VariableNames', {'ShockFrom','FEVD_Share'});
disp(fevdTable)
fprintf('(Note: generalized FEVD does not necessarily sum to 1 across shocks\n');
fprintf('unless residual covariance is diagonal — this is expected, not an error.)\n');

%% ---- GENERALIZED IRF (Pesaran-Shin) ----
irfResp = irf(EstMdl, 'NumObs', 12, 'Method', 'generalized');

save('india_c1_results.mat', 'EstMdl', 'Residuals', 'D', 'irfResp', 'Summary_gc', 'lmarTable');
fprintf('\nSaved india_c1_results.mat\n');


%% ================================================================
%  LOCAL FUNCTION: multivariate LM test for residual autocorrelation
%  Exact formula from Stata's varlmar documentation (Johansen 1995):
%    LM_s = (T - d - 0.5) * ln(|Sigma_hat| / |Sigma_tilde_s|)
%    d = number of coefficients in the AUGMENTED model at lag s
%    df = K^2 (fixed across all tested lags)
%  ================================================================
function tbl = local_varlmar(Y, p, U, mlagmax)
[T, K] = size(Y);
T_eff = T - p;

X0 = ones(T_eff, 1);
for j = 1:p
    X0 = [X0, Y(p+1-j:T-j, :)];
end

Sigma_hat = (U' * U) / T_eff;
detSigma_hat = det(Sigma_hat);

LM = nan(mlagmax,1);
pval = nan(mlagmax,1);
dfout = nan(mlagmax,1);

for s = 1:mlagmax
    % SINGLE lag-s residual term only (not cumulative 1..s)
    Us = zeros(T_eff, K);
    if s < T_eff
        Us(s+1:end, :) = U(1:end-s, :);
    end

    Z = [X0, Us];
    Beta = Z \ Y(p+1:T, :);
    E_s = Y(p+1:T, :) - Z * Beta;
    Sigma_tilde = (E_s' * E_s) / T_eff;

    d = K * size(Z, 2);   % total coefficients across all K equations

    LM(s) = (T_eff - d - 0.5) * log(detSigma_hat / det(Sigma_tilde));
    dfout(s) = K^2;
    pval(s) = 1 - chi2cdf(LM(s), K^2);
end

tbl = table((1:mlagmax)', LM, dfout, pval, ...
    'VariableNames', {'Lag','LM_stat','df','pValue'});
end