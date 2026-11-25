%% ================================================================
%  CORRECTION — India C4+T VECM, rank corrected from 2 to 1
%  The ECT stationarity check showed the 2nd cointegrating relation
%  was spurious (ADF p=0.118). True rank is 1, not 2.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

name = 'India C4+T (M3) [CORRECTED, rank=1]';
allVars = {'ln_nifty', 'ln_india_ip', 'ln_india_m3', 'ln_msci_world', 'us3m_tbills'};
shortRunLags = 2;
rank = 1;   % corrected from 2
respVar = 'ln_india_ip';

fprintf('%s | VECM, rank=%d, short-run lags=%d\n', name, rank, shortRunLags);

Y = data{:, allVars};
validRows = all(~isnan(Y), 2);
Y = Y(validRows, :);
fprintf('N = %d\n\n', size(Y,1));

Mdl = vecm(5, rank, shortRunLags);
Mdl.SeriesNames = allVars;
EstMdl = estimate(Mdl, Y);
summarize(EstMdl)

%% ---- Re-verify: check the (single) ECT is now stationary ----
B = EstMdl.Cointegration;
c0 = EstMdl.CointegrationConstant';
ECT = Y * B + c0;
[h, pval, stat] = adftest(ECT, 'lags', 4, 'model', 'ARD');
if h == 1
    ectVerdict = "Mean-reverting, cointegration confirmed";
else
    ectVerdict = "STILL fails -- needs further investigation";
end
fprintf('\nCorrected ECT check: ADF stat=%.4f, p=%.4f -> %s\n', stat, pval, ectVerdict);

%% ---- Corrected FEVD ----
fprintf('\n=== CORRECTED Generalized FEVD, response = %s, horizon 12 ===\n', respVar);
D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');
respIdx = find(strcmp(allVars, respVar));
fevd_month12 = squeeze(D(12, :, respIdx));
fevdTable = table(allVars', fevd_month12', 'VariableNames', {'ShockFrom','FEVD_Share'});
disp(fevdTable)

fprintf('\nFor comparison, the (incorrect) rank=2 result was:\n');
fprintf('Nifty 31.2%%, IP(own) 57.8%%, M3 0.3%%, MSCI 5.9%%, T-bill 8.6%%\n');

%% ---- Save corrected results ----
save('vecm_India_C4_T_M3_CORRECTED.mat', 'EstMdl', 'D');
fprintf('\nSaved vecm_India_C4_T_M3_CORRECTED.mat\n');