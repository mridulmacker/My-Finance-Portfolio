%% ================================================================
%  STRUCTURAL BREAK TESTING — Chow test at 2008 GFC and COVID-19
%  Tests India's and Korea's stock-return equations (matches the
%  original Stata attempt) for a break in regression coefficients.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

%% ---- Build the equation: stock return on 2 lags of all 4 primary-spec vars ----
% India: dln_nifty on L1/L2 of dln_nifty, dln_india_ip, dln_usd_inr, dln_msci_world
buildEquationTable = @(tbl, stockVar, ipVar, domVar, globVar) ...
    table(tbl.date(3:end), tbl.(stockVar)(3:end), ...
    tbl.(stockVar)(2:end-1), tbl.(stockVar)(1:end-2), ...
    tbl.(ipVar)(2:end-1),    tbl.(ipVar)(1:end-2), ...
    tbl.(domVar)(2:end-1),   tbl.(domVar)(1:end-2), ...
    tbl.(globVar)(2:end-1),  tbl.(globVar)(1:end-2), ...
    'VariableNames', {'date','y','y_L1','y_L2','ip_L1','ip_L2', ...
    'dom_L1','dom_L2','glob_L1','glob_L2'});

%% ================= INDIA =================
Tind = buildEquationTable(data, 'dln_nifty','dln_india_ip','dln_usd_inr','dln_msci_world');
Tind = Tind(~any(ismissing(Tind),2), :);   % drop any remaining NaN rows

fprintf('=== INDIA: stock-return equation, structural break tests ===\n');
fprintf('Sample: %s to %s (%d obs)\n\n', datestr(Tind.date(1)), datestr(Tind.date(end)), height(Tind));

% --- 2008 GFC break (September 2008) ---
bp_gfc = find(Tind.date >= datetime(2008,9,1), 1);
fprintf('--- 2008 GFC break test (breakpoint = %s) ---\n', datestr(Tind.date(bp_gfc)));
chowtest(Tind, bp_gfc, 'Display', 'summary', ...
    'ResponseVariable', 'y', ...
    'PredictorVariables', {'y_L1','y_L2','ip_L1','ip_L2','dom_L1','dom_L2','glob_L1','glob_L2'});

% --- COVID-19 break (March 2020) ---
bp_covid = find(Tind.date >= datetime(2020,3,1), 1);
fprintf('\n--- COVID-19 break test (breakpoint = %s) ---\n', datestr(Tind.date(bp_covid)));
chowtest(Tind, bp_covid, 'Display', 'summary', ...
    'ResponseVariable', 'y', ...
    'PredictorVariables', {'y_L1','y_L2','ip_L1','ip_L2','dom_L1','dom_L2','glob_L1','glob_L2'});

% --- CUSUM test (no break date assumed — detects instability anywhere) ---
fprintf('\n--- CUSUM test (India, no break date assumed) ---\n');
cusumtest(Tind, 'ResponseVariable', 'y', ...
    'PredictorVariables', {'y_L1','y_L2','ip_L1','ip_L2','dom_L1','dom_L2','glob_L1','glob_L2'}, ...
    'Display', 'summary');


%% ================= KOREA =================
Tkor = buildEquationTable(data, 'dln_kospi','dln_korea_ip','dln_usd_kwr','dln_msci_world');
Tkor = Tkor(~any(ismissing(Tkor),2), :);

fprintf('\n\n=== KOREA: stock-return equation, structural break tests ===\n');
fprintf('Sample: %s to %s (%d obs)\n\n', datestr(Tkor.date(1)), datestr(Tkor.date(end)), height(Tkor));

bp_gfc_k = find(Tkor.date >= datetime(2008,9,1), 1);
fprintf('--- 2008 GFC break test (breakpoint = %s) ---\n', datestr(Tkor.date(bp_gfc_k)));
chowtest(Tkor, bp_gfc_k, 'Display', 'summary', ...
    'ResponseVariable', 'y', ...
    'PredictorVariables', {'y_L1','y_L2','ip_L1','ip_L2','dom_L1','dom_L2','glob_L1','glob_L2'});

bp_covid_k = find(Tkor.date >= datetime(2020,3,1), 1);
fprintf('\n--- COVID-19 break test (breakpoint = %s) ---\n', datestr(Tkor.date(bp_covid_k)));
chowtest(Tkor, bp_covid_k, 'Display', 'summary', ...
    'ResponseVariable', 'y', ...
    'PredictorVariables', {'y_L1','y_L2','ip_L1','ip_L2','dom_L1','dom_L2','glob_L1','glob_L2'});

fprintf('\n--- CUSUM test (Korea, no break date assumed) ---\n');
cusumtest(Tkor, 'ResponseVariable', 'y', ...
    'PredictorVariables', {'y_L1','y_L2','ip_L1','ip_L2','dom_L1','dom_L2','glob_L1','glob_L2'}, ...
    'Display', 'summary');