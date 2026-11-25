%% ================================================================
%  ERROR-CORRECTION TERM CHECK (per Ben's feedback, point 1d)
%  Extracts the cointegrating relation(s) from each fitted VECM,
%  computes the implied error-correction term, plots it, and runs
%  an ADF test to formally check it's actually stationary/mean-reverting.
%  Requires the vecm_*.mat files already saved from step27.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

specs = {
    'India C1+T (FX)',  'ln_nifty', 'ln_india_ip', 'ln_usd_inr',     'ln_msci_world', 'us3m_tbills';
    'India C2+T (CPI)', 'ln_nifty', 'ln_india_ip', 'ln_india_cpi',   'ln_msci_world', 'us3m_tbills';
    'India C3+T (IR)',  'ln_nifty', 'ln_india_ip', 'india_call_rate','ln_msci_world', 'us3m_tbills';
    'India C4+T (M3)',  'ln_nifty', 'ln_india_ip', 'ln_india_m3',    'ln_msci_world', 'us3m_tbills';
    'Korea D1+T (FX)',  'ln_kospi', 'ln_korea_ip', 'ln_usd_kwr',     'ln_msci_world', 'us3m_tbills';
    'Korea D2+T (CPI)', 'ln_kospi', 'ln_korea_ip', 'ln_korea_cpi',   'ln_msci_world', 'us3m_tbills';
    'Korea D3+T (IR)',  'ln_kospi', 'ln_korea_ip', 'korea_ir',       'ln_msci_world', 'us3m_tbills';
    'Korea D4+T (M3)',  'ln_kospi', 'ln_korea_ip', 'ln_korea_m3',    'ln_msci_world', 'us3m_tbills';
    };

fprintf('=== ERROR-CORRECTION TERM STATIONARITY CHECK ===\n');
fprintf('If the ECT fails its own ADF test, the cointegrating relation is\n');
fprintf('likely spurious regardless of what the trace/maxeig tests showed.\n\n');

ectSummary = {};

for i = 1:size(specs,1)
    name = specs{i,1};
    allVars = specs(i,2:6);

    cleanName = regexprep(name, '[^a-zA-Z0-9]', '_');
    vecmFile = sprintf('vecm_%s.mat', cleanName);

    fprintf('\n----------------------------------------------------\n');
    fprintf('%s\n', name);
    fprintf('----------------------------------------------------\n');

    try
        S = load(vecmFile);
        EstMdl = S.EstMdl;
    catch ME
        fprintf('Could not load %s: %s\n', vecmFile, ME.message);
        continue;
    end

    Y = data{:, allVars};
    validRows = all(~isnan(Y), 2);
    Y = Y(validRows, :);
    theseDates = data.date(validRows);

    B = EstMdl.Cointegration;   % n x rank matrix
    rank = size(B, 2);

    if isprop(EstMdl, 'CointegrationConstant') && ~isempty(EstMdl.CointegrationConstant)
        c0 = EstMdl.CointegrationConstant';   % 1 x rank
    else
        c0 = zeros(1, rank);
    end

    % Error-correction term(s): ECT_t = Y_t * B + c0, for each cointegrating relation
    ECT = Y * B + c0;   % T x rank

    figure('Name', name);
    for r = 1:rank
        subplot(rank, 1, r);
        plot(theseDates, ECT(:,r), 'LineWidth', 1.2);
        yline(0, '--k');
        title(sprintf('%s -- Error-Correction Term %d of %d', name, r, rank));
        xlabel('Date'); ylabel('ECT value');
        grid on;
    end

    % Formal ADF test on each ECT series (should be stationary if genuine)
    for r = 1:rank
        ectSeries = ECT(:,r);
        [h, pval, stat] = adftest(ectSeries, 'lags', 4, 'model', 'ARD');
        verdict = "NOT mean-reverting (fails ADF) -- cointegration doubtful";
        if h == 1
            verdict = "Mean-reverting (passes ADF) -- supports cointegration";
        end
        fprintf('ECT %d/%d: ADF stat=%.4f, p=%.4f -> %s\n', r, rank, stat, pval, verdict);
        ectSummary(end+1,:) = {name, r, pval, verdict};
    end
end

fprintf('\n\n==================== SUMMARY: ECT stationarity across all specs ====================\n');
summaryTable = cell2table(ectSummary, 'VariableNames', {'Spec','ECT_number','pValue','Verdict'});
disp(summaryTable)