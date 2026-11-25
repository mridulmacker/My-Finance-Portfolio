%% ================================================================
%  JOHANSEN COINTEGRATION TEST — 8 T-bill specs, 5-variable systems
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

% Spec table: {name, stockLevelVar, ipLevelVar, domesticLevelVar, globalLevelVar, tbillLevelVar, lag}
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

fprintf('=== JOHANSEN COINTEGRATION TESTS (5-variable systems, with T-bill, levels) ===\n');
fprintf('Model: H1 (constant in cointegrating relation, no trend)\n\n');

resultsSummary = {};

for i = 1:size(specs,1)
    name = specs{i,1};
    v1 = specs{i,2}; v2 = specs{i,3}; v3 = specs{i,4}; v4 = specs{i,5}; v5 = specs{i,6};
    lag = specs{i,7};

    Y = data{:, {v1,v2,v3,v4,v5}};
    validRows = all(~isnan(Y), 2);
    Y = Y(validRows, :);

    fprintf('\n--------------------------------------------------\n');
    fprintf('%s | Variables: %s, %s, %s, %s, %s | Lags: %d | N: %d\n', ...
        name, v1, v2, v3, v4, v5, lag, size(Y,1));
    fprintf('--------------------------------------------------\n');

    [h, pValue, stat, cValue] = jcitest(Y, 'Model', 'H1', 'Lags', lag, ...
        'Test', {'trace'}, 'Display', 'summary');

    % Determine rank via sequential procedure: first non-rejection stops the count
    rank = 0;
    fieldNames = {'r0','r1','r2','r3','r4'};
    for rIdx = 1:4  % max possible rank for a 5-variable system is 4
        if h.(fieldNames{rIdx})(1) == true
            rank = rIdx;  % reject "rank <= rIdx-1", so at least rIdx cointegrating relations
        else
            break;  % first non-rejection stops the sequential procedure
        end
    end

    resultsSummary(end+1,:) = {name, rank, pValue.r0(1)};
end

fprintf('\n\n==================== SUMMARY (all 8 T-bill specs) ====================\n');
summaryTable = cell2table(resultsSummary, 'VariableNames', {'Spec','Rank_chosen','pValue_r0'});
disp(summaryTable)