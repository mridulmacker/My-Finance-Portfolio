%% ================================================================
%  TABLE C2 — RESIDUAL DIAGNOSTICS, all 8 five-variable specs
%  Prints one summary row per specification (system-level where possible,
%  worst-equation where the test is per-equation).
% ================================================================
clear; clc;
load('master_dataset_transformed.mat');

specs = {
    'India FX',  'dln_nifty','dln_india_ip','dln_usd_inr',     'dln_msci_world','d_us3mtbills', 2;
    'India CPI', 'dln_nifty','dln_india_ip','dln_india_cpi',   'dln_msci_world','d_us3mtbills', 2;
    'India IR',  'dln_nifty','dln_india_ip','d_indiacallrate', 'dln_msci_world','d_us3mtbills', 5;
    'India M3',  'dln_nifty','dln_india_ip','dln_india_m3',    'dln_msci_world','d_us3mtbills', 2;
    'Korea FX',  'dln_kospi','dln_korea_ip','dln_usd_kwr',     'dln_msci_world','d_us3mtbills', 2;
    'Korea CPI', 'dln_kospi','dln_korea_ip','dln_korea_cpi',   'dln_msci_world','d_us3mtbills', 2;
    'Korea IR',  'dln_kospi','dln_korea_ip','d_koreair',       'dln_msci_world','d_us3mtbills', 2;
    'Korea M3',  'dln_kospi','dln_korea_ip','dln_korea_m3',    'dln_msci_world','d_us3mtbills', 3;
    };

fprintf('SPEC | Stable(maxEig) | LM_min_p(12) | JB_all_reject | ARCH_min_p(4)\n');
fprintf('-------------------------------------------------------------------\n');
for i = 1:size(specs,1)
    name = specs{i,1}; vars = specs(i,2:6); p = specs{i,7};
    Y = data{:, vars}; Y = Y(all(~isnan(Y),2), :);
    Mdl = varm(5, p);
    EstMdl = estimate(Mdl, Y);

    % Stability: largest eigenvalue modulus of companion matrix
    % (summarize gives it, but compute directly for safety)
    ar = EstMdl.AR; k = numel(ar); n = 5;
    C = zeros(n*k);
    for L=1:k, C(1:n,(L-1)*n+1:L*n) = ar{L}; end
    if k>1, C(n+1:end,1:n*(k-1)) = eye(n*(k-1)); end
    maxEig = max(abs(eig(C)));

    E = infer(EstMdl, Y);   % residuals

    % Per-equation: LM (Ljung-Box 12), JB, ARCH(4). Summarise across equations.
    lb_ps = zeros(1,5); jb_rej = 0; arch_ps = zeros(1,5);
    for q = 1:5
        e = E(:,q);
        [~, lb_ps(q)] = lbqtest(e, 'Lags', 12);
        [hj, ~] = jbtest(e); jb_rej = jb_rej + hj;   % count equations rejecting normality
        [~, arch_ps(q)] = archtest(e, 'Lags', 4);
    end
    fprintf('%-9s | %s (%.3f) | %.3f | %d/5 | %.3f\n', ...
        name, string(maxEig<1), maxEig, min(lb_ps), jb_rej, min(arch_ps));
end
fprintf('\nReading: Stable=true if maxEig<1. LM_min_p = smallest Ljung-Box p across the 5 equations\n');
fprintf('(low = autocorrelation in at least one equation). JB_all_reject = how many of 5 equations\n');
fprintf('reject normality. ARCH_min_p = smallest ARCH p across equations (low = ARCH present).\n');