%% ================================================================
%  INDIA C1 — REVERSE-DIRECTION FEVD (response = Nifty, not IP)
%  Completes the symmetry check already done for Korea D1.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

vars = {'dln_nifty','dln_india_ip','dln_usd_inr','dln_msci_world'};
Y = data{:, vars};
validRows = all(~isnan(Y), 2);
Y = Y(validRows, :);

p = 2;  % validated lag for India C1
Mdl = varm(4, p);
Mdl.SeriesNames = vars;
[EstMdl, ~, ~, ~] = estimate(Mdl, Y);

D = fevd(EstMdl, 'NumObs', 12, 'Method', 'generalized');

fprintf('=== Generalized FEVD, response = Nifty (reverse direction), horizon 12 ===\n');
respIdx = find(strcmp(vars, 'dln_nifty'));
fevd_month12 = squeeze(D(12, :, respIdx));
fevdTable = table(vars', fevd_month12', 'VariableNames', {'ShockFrom','FEVD_Share'});
disp(fevdTable)

fprintf('\nFor reference, forward direction (response = India IP):\n');
respIdx2 = find(strcmp(vars, 'dln_india_ip'));
fevd_month12_ip = squeeze(D(12, :, respIdx2));
fevdTable2 = table(vars', fevd_month12_ip', 'VariableNames', {'ShockFrom','FEVD_Share'});
disp(fevdTable2)