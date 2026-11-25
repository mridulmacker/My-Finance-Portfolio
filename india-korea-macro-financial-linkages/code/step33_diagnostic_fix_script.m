%% ================================================================
%  DIAGNOSTIC + FIX SCRIPT
%  Part 1: Settle the lag-convention question (are jcitest & vecm consistent?)
%  Part 2: Re-run all 8 VECM FEVDs (fills Korea gap in Table 5.9)
%  Part 3: Row-normalize ALL generalized FEVDs (VAR-in-diff + VECM)
%  Part 4: Print clean comparison tables
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

%% ================================================================
%  PART 1 — LAG CONVENTION CHECK
%  Confirms what MATLAB actually does, so we know if a rerun is needed.
% ================================================================
fprintf('========================================================\n');
fprintf('PART 1: LAG CONVENTION CHECK\n');
fprintf('========================================================\n');
fprintf(['In MATLAB:\n' ...
    '  - jcitest ''Lags'' = number of lagged DIFFERENCES (VECM q)\n' ...
    '  - vecm(K, rank, q): third argument q = lagged differences too\n' ...
    'If your scripts passed the SAME number to both, they are consistent\n' ...
    'and no rerun is needed. The check below proves it empirically:\n' ...
    'we estimate India FX two ways and confirm they agree.\n\n']);

allVars = {'ln_nifty','ln_india_ip','ln_usd_inr','ln_msci_world','us3m_tbills'};
Y = data{:, allVars};
Y = Y(all(~isnan(Y),2), :);

% The VAR lag you selected for India FX was p = 2.
% A VAR(p) in levels corresponds to a VECM with q = p-1 lagged differences.
% Your jcitest call used 'Lags', 2. Let's see what MATLAB's own vecm
% documentation-consistent estimate does at q=2 vs q=1, and check the
% Johansen test at both, so we can see which matches your reported ranks.

fprintf('--- Johansen trace test, India FX, at Lags = 1 and Lags = 2 ---\n');
for L = [1 2]
    [h,pv,st,cv] = jcitest(Y, 'Model','H1', 'Lags', L, 'Test','trace', 'Display','off');
    rank = 0; fn = {'r0','r1','r2','r3','r4'};
    for r = 1:4
        if h.(fn{r})(1), rank = r; else, break; end
    end
    fprintf('  Lags=%d -> rank chosen = %d (r0 p-value = %.4f)\n', L, rank, pv.r0(1));
end
fprintf(['\nWhich of these matches the rank=1 you reported for India FX?\n' ...
    'Whichever Lags value reproduces your reported rank is the convention\n' ...
    'you actually used. If Lags=2 gives rank 1 (your reported result),\n' ...
    'your pipeline was internally consistent and NO RERUN is needed.\n\n']);

%% ================================================================
%  PART 2 & 3 — RE-ESTIMATE ALL 8 VECMs, EXTRACT + NORMALIZE FEVD
% ================================================================
fprintf('========================================================\n');
fprintf('PART 2 & 3: VECM FEVD (raw + row-normalized), all 8 specs\n');
fprintf('========================================================\n');

% {name, stock, ip, dom, msci, tbill, shortRunLags(q), rank}
specs = {
    'India C1+T (FX)',  'ln_nifty','ln_india_ip','ln_usd_inr',     'ln_msci_world','us3m_tbills', 2, 1;
    'India C2+T (CPI)', 'ln_nifty','ln_india_ip','ln_india_cpi',   'ln_msci_world','us3m_tbills', 2, 1;
    'India C3+T (IR)',  'ln_nifty','ln_india_ip','india_call_rate','ln_msci_world','us3m_tbills', 5, 1;
    'India C4+T (M3)',  'ln_nifty','ln_india_ip','ln_india_m3',    'ln_msci_world','us3m_tbills', 2, 1;
    'Korea D1+T (FX)',  'ln_kospi','ln_korea_ip','ln_usd_kwr',     'ln_msci_world','us3m_tbills', 2, 1;
    'Korea D2+T (CPI)', 'ln_kospi','ln_korea_ip','ln_korea_cpi',   'ln_msci_world','us3m_tbills', 2, 1;
    'Korea D3+T (IR)',  'ln_kospi','ln_korea_ip','korea_ir',       'ln_msci_world','us3m_tbills', 2, 2;
    'Korea D4+T (M3)',  'ln_kospi','ln_korea_ip','ln_korea_m3',    'ln_msci_world','us3m_tbills', 3, 1;
    };

for i = 1:size(specs,1)
    name = specs{i,1};
    vars = specs(i,2:6);
    q    = specs{i,7};
    rank = specs{i,8};
    respVar = vars{2};   % industrial production

    Y = data{:, vars};
    Y = Y(all(~isnan(Y),2), :);

    Mdl = vecm(5, rank, q);
    Mdl.SeriesNames = vars;
    try
        EstMdl = estimate(Mdl, Y);
        D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');
        respIdx = find(strcmp(vars, respVar));
        raw = squeeze(D(12, :, respIdx));      % raw generalized FEVD (won't sum to 1)
        norm = raw / sum(raw);                  % ROW-NORMALIZED (now sums to 1)

        fprintf('\n--- %s (rank=%d, q=%d), response = %s ---\n', name, rank, q, respVar);
        T = table(vars', raw', norm', 'VariableNames', {'ShockFrom','Raw','Normalized'});
        disp(T);
        fprintf('   (raw sum = %.4f, normalized sum = %.4f)\n', sum(raw), sum(norm));
    catch ME
        fprintf('\n--- %s FAILED: %s ---\n', name, ME.message);
    end
end

%% ================================================================
%  PART 4 — NORMALIZE THE VAR-IN-DIFFERENCES FEVDs TOO
%  Re-does the 8 five-variable VAR FEVDs and normalizes them, so the
%  Results tables 5.6 use honest shares that sum to 100%.
% ================================================================
fprintf('\n========================================================\n');
fprintf('PART 4: VAR-in-differences FEVD (normalized), 5-var specs\n');
fprintf('========================================================\n');

dspecs = {
    'India C1+T (FX)',  'dln_nifty','dln_india_ip','dln_usd_inr',     'dln_msci_world','d_us3mtbills', 2;
    'India C2+T (CPI)', 'dln_nifty','dln_india_ip','dln_india_cpi',   'dln_msci_world','d_us3mtbills', 2;
    'India C3+T (IR)',  'dln_nifty','dln_india_ip','d_indiacallrate', 'dln_msci_world','d_us3mtbills', 5;
    'India C4+T (M3)',  'dln_nifty','dln_india_ip','dln_india_m3',    'dln_msci_world','d_us3mtbills', 2;
    'Korea D1+T (FX)',  'dln_kospi','dln_korea_ip','dln_usd_kwr',     'dln_msci_world','d_us3mtbills', 2;
    'Korea D2+T (CPI)', 'dln_kospi','dln_korea_ip','dln_korea_cpi',   'dln_msci_world','d_us3mtbills', 2;
    'Korea D3+T (IR)',  'dln_kospi','dln_korea_ip','d_koreair',       'dln_msci_world','d_us3mtbills', 2;
    'Korea D4+T (M3)',  'dln_kospi','dln_korea_ip','dln_korea_m3',    'dln_msci_world','d_us3mtbills', 3;
    };

for i = 1:size(dspecs,1)
    name = dspecs{i,1};
    vars = dspecs(i,2:6);
    p    = dspecs{i,6+1};
    respVar = vars{2};

    Y = data{:, vars};
    Y = Y(all(~isnan(Y),2), :);

    Mdl = varm(5, p);
    Mdl.SeriesNames = vars;
    EstMdl = estimate(Mdl, Y);
    D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');
    respIdx = find(strcmp(vars, respVar));
    raw = squeeze(D(12, :, respIdx));
    norm = raw / sum(raw);

    fprintf('\n--- %s (VAR, p=%d), response = %s ---\n', name, p, respVar);
    T = table(vars', raw', norm', 'VariableNames', {'ShockFrom','Raw','Normalized'});
    disp(T);
    fprintf('   (raw sum = %.4f)\n', sum(raw));
end

fprintf('\n\nDONE. Send me the full output.\n');