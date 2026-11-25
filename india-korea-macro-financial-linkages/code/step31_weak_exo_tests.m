%% ================================================================
%  WEAK EXOGENEITY TESTS — CORRECTED syntax, using ACon restriction matrix
%  H0: R'*alpha = 0, where R selects the row(s) of alpha being restricted.
%  Joint test (MSCI+Tbill): R has 2 columns, each a selector for one row.
%  Stock test: R has 1 column, selector for the stock row.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

specs = {
    'India C1+T (FX)',  'ln_nifty', 'ln_india_ip', 'ln_usd_inr',     'ln_msci_world', 'us3m_tbills', 2, 1, '';
    'India C2+T (CPI)', 'ln_nifty', 'ln_india_ip', 'ln_india_cpi',   'ln_msci_world', 'us3m_tbills', 2, 1, '';
    'India C3+T (IR)',  'ln_nifty', 'ln_india_ip', 'india_call_rate','ln_msci_world', 'us3m_tbills', 5, 1, 'CAVEAT: rank not clearly supported - illustrative only';
    'India C4+T (M3)',  'ln_nifty', 'ln_india_ip', 'ln_india_m3',    'ln_msci_world', 'us3m_tbills', 2, 1, 'Corrected rank (was 2)';
    'Korea D1+T (FX)',  'ln_kospi', 'ln_korea_ip', 'ln_usd_kwr',     'ln_msci_world', 'us3m_tbills', 2, 1, '';
    'Korea D2+T (CPI)', 'ln_kospi', 'ln_korea_ip', 'ln_korea_cpi',   'ln_msci_world', 'us3m_tbills', 2, 1, 'CAVEAT: 2/5 deterministic forms found no cointegration';
    'Korea D3+T (IR)',  'ln_kospi', 'ln_korea_ip', 'korea_ir',       'ln_msci_world', 'us3m_tbills', 2, 2, 'CAVEAT: small-sample correction borderline';
    'Korea D4+T (M3)',  'ln_kospi', 'ln_korea_ip', 'ln_korea_m3',    'ln_msci_world', 'us3m_tbills', 3, 1, '';
    };

results = {};

for i = 1:size(specs,1)
    name = specs{i,1};
    allVars = specs(i,2:6);
    lag = specs{i,7};
    rank = specs{i,8};
    caveat = specs{i,9};

    fprintf('\n====================================================\n');
    fprintf('%s | rank=%d, lag=%d\n', name, rank, lag);
    if ~isempty(caveat), fprintf('%s\n', caveat); end
    fprintf('====================================================\n');

    Y = data{:, allVars};
    validRows = all(~isnan(Y), 2);
    Y = Y(validRows, :);
    n = numel(allVars);   % = 5 always here: [stock, IP, domestic, MSCI, Tbill]

    % ---- Test 1: JOINT — MSCI World (row 4) + T-bill (row 5) both = 0 ----
    R_joint = zeros(n, 2);
    R_joint(4,1) = 1;   % selects MSCI row
    R_joint(5,2) = 1;   % selects T-bill row
    try
        [h_g, p_g, stat_g, cV_g] = jcontest(Y, rank, "ACon", R_joint, ...
            'Model', 'H1', 'Lags', lag);
        fprintf('\n--- Global factors (MSCI+Tbill) jointly weakly exogenous? ---\n');
        fprintf('h=%d, p=%.4f, LR stat=%.4f, cValue=%.4f\n', h_g, p_g, stat_g, cV_g);
    catch ME
        fprintf('\nJOINT TEST FAILED: %s\n', ME.message);
        h_g = NaN; p_g = NaN; stat_g = NaN;
    end

    % ---- Test 2: domestic stock index (row 1) = 0 ----
    R_stock = zeros(n, 1);
    R_stock(1,1) = 1;   % selects stock row
    try
        [h_s, p_s, stat_s, cV_s] = jcontest(Y, rank, "ACon", R_stock, ...
            'Model', 'H1', 'Lags', lag);
        fprintf('\n--- Domestic stock index weakly exogenous? ---\n');
        fprintf('h=%d, p=%.4f, LR stat=%.4f, cValue=%.4f\n', h_s, p_s, stat_s, cV_s);
    catch ME
        fprintf('\nSTOCK TEST FAILED: %s\n', ME.message);
        h_s = NaN; p_s = NaN; stat_s = NaN;
    end

    results(end+1,:) = {name, rank, h_g, p_g, stat_g, h_s, p_s, stat_s, caveat};
end

fprintf('\n\n==================== SUMMARY: WEAK EXOGENEITY (all 8 specs) ====================\n');
summaryTable = cell2table(results, 'VariableNames', ...
    {'Spec','Rank','Global_h','Global_p','Global_LR','Stock_h','Stock_p','Stock_LR','Caveat'});
disp(summaryTable)

% h = 1 -> REJECT weak exogeneity (variable DOES adjust -- not weakly exogenous)
% h = 0 -> FAIL TO REJECT (consistent with weak exogeneity -- variable "drives")