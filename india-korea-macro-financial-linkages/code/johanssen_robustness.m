%% ================================================================
%  JOHANSEN ROBUSTNESS CHECKS (per Ben's feedback)
%  1. Max-eigenvalue test alongside trace test
%  2. Sensitivity to deterministic assumption (all 5 Johansen forms)
%  3. Small-sample (Reimers 1992) correction to the trace statistic
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

models = {'H2','H1*','H1','H*','H'};

for i = 1:size(specs,1)
    name = specs{i,1};
    allVars = specs(i,2:6);
    lag = specs{i,7};

    fprintf('\n\n====================================================\n');
    fprintf('%s\n', name);
    fprintf('====================================================\n');

    Y = data{:, allVars};
    validRows = all(~isnan(Y), 2);
    Y = Y(validRows, :);
    T = size(Y,1) - lag;   % effective sample size
    n = numel(allVars);    % number of variables

    %% ---- 1 & 2: trace AND maxeig, across all 5 deterministic forms ----
    fprintf('\n--- Rank decision (r=0 test) across deterministic forms and test types ---\n');
    fprintf('%-8s %-10s %-10s %-10s\n', 'Model', 'Trace_p', 'MaxEig_p', 'Rank(trace)');

    for m = 1:numel(models)
        try
            [h_tr, p_tr, stat_tr, cV_tr] = jcitest(Y, 'Model', models{m}, 'Lags', lag, ...
                'Test', {'trace'}, 'Display', 'off');
            [h_me, p_me, stat_me, cV_me] = jcitest(Y, 'Model', models{m}, 'Lags', lag, ...
                'Test', {'maxeig'}, 'Display', 'off');

            % Determine rank via sequential procedure on trace test
            rank = 0;
            fieldNames = {'r0','r1','r2','r3','r4'};
            for rIdx = 1:(n-1)
                if h_tr.(fieldNames{rIdx})(1) == true
                    rank = rIdx;
                else
                    break;
                end
            end

            fprintf('%-8s %-10.4f %-10.4f %-10d\n', models{m}, p_tr.r0(1), p_me.r0(1), rank);
        catch ME
            fprintf('%-8s FAILED: %s\n', models{m}, ME.message);
        end
    end

    %% ---- 3: Reimers (1992) small-sample correction, on the H1 baseline ----
    fprintf('\n--- Small-sample (Reimers 1992) correction, Model=H1 ---\n');
    try
        [h_base, p_base, stat_base, cV_base] = jcitest(Y, 'Model', 'H1', 'Lags', lag, ...
            'Test', {'trace'}, 'Display', 'off');

        correction_factor = (T - n*lag) / T;
        fprintf('Effective N=%d, variables=%d, lags=%d, correction factor=%.4f\n', T, n, lag, correction_factor);
        fprintf('%-6s %-12s %-12s %-12s %-10s\n', 'rank', 'raw_stat', 'corrected', 'cValue', 'still_reject?');

        fieldNames = {'r0','r1','r2','r3','r4'};
        for rIdx = 1:n
            raw = stat_base.(fieldNames{rIdx})(1);
            corrected = raw * correction_factor;
            cv = cV_base.(fieldNames{rIdx})(1);
            stillReject = corrected > cv;
            fprintf('%-6d %-12.4f %-12.4f %-12.4f %-10s\n', rIdx-1, raw, corrected, cv, string(stillReject));
        end
    catch ME
        fprintf('Correction step failed: %s\n', ME.message);
    end
end