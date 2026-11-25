%% ================================================================
%  VECM ESTIMATION — all 8 T-bill (5-variable) specs
%  Ranks confirmed via Johansen cointegration test (step26)
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

% Spec table: {name, stockVar, ipVar, domVar, msciVar, tbillVar, shortRunLags, rank}
specs = {
    'India C1+T (FX)',  'ln_nifty', 'ln_india_ip', 'ln_usd_inr',     'ln_msci_world', 'us3m_tbills', 2, 1;
    'India C2+T (CPI)', 'ln_nifty', 'ln_india_ip', 'ln_india_cpi',   'ln_msci_world', 'us3m_tbills', 2, 1;
    'India C3+T (IR)',  'ln_nifty', 'ln_india_ip', 'india_call_rate','ln_msci_world', 'us3m_tbills', 5, 1;
    'India C4+T (M3)',  'ln_nifty', 'ln_india_ip', 'ln_india_m3',    'ln_msci_world', 'us3m_tbills', 2, 2;
    'Korea D1+T (FX)',  'ln_kospi', 'ln_korea_ip', 'ln_usd_kwr',     'ln_msci_world', 'us3m_tbills', 2, 1;
    'Korea D2+T (CPI)', 'ln_kospi', 'ln_korea_ip', 'ln_korea_cpi',   'ln_msci_world', 'us3m_tbills', 2, 1;
    'Korea D3+T (IR)',  'ln_kospi', 'ln_korea_ip', 'korea_ir',       'ln_msci_world', 'us3m_tbills', 2, 2;
    'Korea D4+T (M3)',  'ln_kospi', 'ln_korea_ip', 'ln_korea_m3',    'ln_msci_world', 'us3m_tbills', 3, 1;
    };

for i = 1:size(specs,1)
    name = specs{i,1};
    v1 = specs{i,2}; v2 = specs{i,3}; v3 = specs{i,4}; v4 = specs{i,5}; v5 = specs{i,6};
    shortRunLags = specs{i,7};
    rank = specs{i,8};
    respVar = v2;
    allVars = {v1,v2,v3,v4,v5};

    fprintf('\n\n====================================================\n');
    fprintf('%s | VECM, rank=%d, short-run lags=%d\n', name, rank, shortRunLags);
    fprintf('====================================================\n');

    Y = data{:, allVars};
    validRows = all(~isnan(Y), 2);
    Y = Y(validRows, :);
    fprintf('N = %d\n\n', size(Y,1));

    Mdl = vecm(5, rank, shortRunLags);
    Mdl.SeriesNames = allVars;
    try
        EstMdl = estimate(Mdl, Y);
    catch ME
        fprintf('ESTIMATION FAILED: %s\n', ME.message);
        continue;
    end
    summarize(EstMdl)

    try
        [h_gc, Summary_gc] = gctest(EstMdl, 'Display', false);
        fprintf('\n=== Granger causality ===\n');
        disp(Summary_gc)
    catch ME
        fprintf('\ngctest on VECM failed/unsupported: %s\n', ME.message);
    end

    try
        D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');
        respIdx = find(strcmp(allVars, respVar));
        fevd_month12 = squeeze(D(12, :, respIdx));
        fevdTable = table(allVars', fevd_month12', 'VariableNames', {'ShockFrom','FEVD_Share'});
        fprintf('\n=== Generalized FEVD, response = %s, horizon 12 ===\n', respVar);
        disp(fevdTable)
    catch ME
        fprintf('\nFEVD on VECM failed: %s\n', ME.message);
    end

    try
        irfResp = irf(EstMdl, 'NumObs', 12, 'Method', 'generalized');
        fprintf('\nIRF computed successfully.\n');
    catch ME
        fprintf('\nIRF on VECM failed: %s\n', ME.message);
    end

    cleanName = regexprep(name, '[^a-zA-Z0-9]', '_');
    saveName = sprintf('vecm_%s.mat', cleanName);
    save(saveName, 'EstMdl');
end