%% ================================================================
%  TODA-YAMAMOTO GRANGER CAUSALITY — all 8 T-bill (5-variable) specs
%  Uses LEVEL data directly (no VECM object needed), robust to
%  cointegration/integration via the Integration name-value argument.
%  Leave-one-out style: tests each variable -> each other variable,
%  conditioning on the remaining 3.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

specs = {
    'India C1+T (FX)',  'ln_nifty', 'ln_india_ip', 'ln_usd_inr',     'ln_msci_world', 'us3m_tbills', 2;
    'India C2+T (CPI)', 'ln_nifty', 'ln_india_ip', 'ln_india_cpi',   'ln_msci_world', 'us3m_tbills', 2;
    'India C3+T (IR)',  'ln_nifty', 'ln_india_ip', 'india_call_rate','ln_msci_world', 'us3m_tbills', 5;
    'India C4+T (M3)',  'ln_nifty', 'ln_india_ip', 'ln_india_m3',    'ln_msci_world', 'us3m_tbills', 2;
    'Korea D1+T (FX)',  'ln_kospi', 'ln_korea_ip', 'ln_usd_kwr',     'ln_msci_world', 'us3m_tbills', 2;
    'Korea D2+T (CPI)', 'ln_kospi', 'ln_korea_ip', 'ln_korea_cpi',   'ln_msci_world', 'us3m_tbills', 2;
    'Korea D3+T (IR)',  'ln_kospi', 'ln_korea_ip', 'korea_ir',       'ln_msci_world', 'us3m_tbills', 2;
    'Korea D4+T (M3)',  'ln_kospi', 'ln_korea_ip', 'ln_korea_m3',    'ln_msci_world', 'us3m_tbills', 3;
    };

for i = 1:size(specs,1)
    name = specs{i,1};
    allVars = specs(i,2:6);
    lag = specs{i,7};

    fprintf('\n\n====================================================\n');
    fprintf('%s | Toda-Yamamoto Granger causality (lags=%d, Integration=1)\n', name, lag);
    fprintf('====================================================\n');

    Y = data{:, allVars};
    validRows = all(~isnan(Y), 2);
    Y = Y(validRows, :);

    n = numel(allVars);
    results = {};

    for src = 1:n
        for tgt = 1:n
            if src == tgt
                continue;
            end
            condIdx = setdiff(1:n, [src, tgt]);
            Y1 = Y(:, src);
            Y2 = Y(:, tgt);
            Y3 = Y(:, condIdx);

            try
                [h, pval] = gctest(Y1, Y2, Y3, 'NumLags', lag, 'Integration', 1, 'Test', 'f');
                verdict = "Cannot reject H0";
                if h == 1
                    verdict = "Reject H0 (causal)";
                end
                results(end+1,:) = {allVars{src}, allVars{tgt}, pval, verdict};
            catch ME
                results(end+1,:) = {allVars{src}, allVars{tgt}, NaN, ['ERROR: ' ME.message]};
            end
        end
    end

    resultsTable = cell2table(results, 'VariableNames', {'Cause','Effect','pValue','Decision'});
    disp(resultsTable)
end