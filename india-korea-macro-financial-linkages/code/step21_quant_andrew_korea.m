%% ================================================================
%  QUANDT-ANDREWS SCAN — Korea stock-return equation
%  Scans every candidate break date (trimming 15% from each end,
%  standard convention) and finds the date of MAXIMUM instability.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

buildEquationTable = @(tbl, stockVar, ipVar, domVar, globVar) ...
    table(tbl.date(3:end), tbl.(stockVar)(3:end), ...
    tbl.(stockVar)(2:end-1), tbl.(stockVar)(1:end-2), ...
    tbl.(ipVar)(2:end-1),    tbl.(ipVar)(1:end-2), ...
    tbl.(domVar)(2:end-1),   tbl.(domVar)(1:end-2), ...
    tbl.(globVar)(2:end-1),  tbl.(globVar)(1:end-2), ...
    'VariableNames', {'date','y','y_L1','y_L2','ip_L1','ip_L2', ...
    'dom_L1','dom_L2','glob_L1','glob_L2'});

Tkor = buildEquationTable(data, 'dln_kospi','dln_korea_ip','dln_usd_kwr','dln_msci_world');
Tkor = Tkor(~any(ismissing(Tkor),2), :);
N = height(Tkor);

% Standard Andrews (1993) trimming: search only the middle 70% of the sample
trimLow  = round(0.15 * N);
trimHigh = round(0.85 * N);

predictors = {'y_L1','y_L2','ip_L1','ip_L2','dom_L1','dom_L2','glob_L1','glob_L2'};

stats = nan(N,1);
pvals = nan(N,1);

fprintf('Scanning breakpoints %d to %d (of %d total obs)...\n', trimLow, trimHigh, N);
for bp = trimLow:trimHigh
    result = chowtest(Tkor, bp, 'ResponseVariable', 'y', ...
        'PredictorVariables', predictors, 'Display', 'off');
    stats(bp) = result.stat;
    pvals(bp) = result.pValue;
end

[maxStat, maxIdx] = max(stats);
fprintf('\n=== Quandt-Andrews result ===\n');
fprintf('Maximum Chow statistic: %.4f at breakpoint %d (%s)\n', maxStat, maxIdx, datestr(Tkor.date(maxIdx)));
fprintf('Point-wise p-value at that date: %.4f\n', pvals(maxIdx));
fprintf(['NOTE: this point p-value is optimistic — it does not account for\n' ...
    'searching over %d candidate dates (data-mining bias). The true\n' ...
    'Quandt-Andrews sup-F test requires Andrews (1993) adjusted critical\n' ...
    'values, which are higher than a single-date Chow critical value.\n' ...
    'Treat this as identifying WHERE instability is concentrated,\n' ...
    'not as a formally significant break on its own.\n'], trimHigh-trimLow+1);

% Plot the statistic over time to visualize where instability concentrates
figure;
plot(Tkor.date(trimLow:trimHigh), stats(trimLow:trimHigh), 'LineWidth', 1.5);
hold on;
yline(maxStat, '--r');
xline(Tkor.date(maxIdx), '--k', datestr(Tkor.date(maxIdx),'mmm yyyy'));
title('Korea: Chow statistic across candidate break dates');
xlabel('Candidate breakpoint date');
ylabel('Chow F-statistic');
grid on;