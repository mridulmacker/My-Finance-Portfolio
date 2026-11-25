%% ================================================================
%  ADDITIONAL DIAGNOSTICS + BOOTSTRAP FEVD
%  Block 1: Residual diagnostics (Jarque-Bera, ARCH, Ljung-Box) — all 8 specs
%  Block 2: Bootstrap CIs on FEVD of IP + test India vs Korea global-share gap
%  Block 3: COVID-dummy robustness on the key India T-bill -> IP result
%  ================================================================
clear; clc; rng(42);   % fixed seed for reproducible bootstrap
load('master_dataset_transformed.mat');

% {name, stock, ip, dom, msci, tbill, p(VARlag)}  -- differenced (model) forms
specs = {
 'India C1+T (FX)',  'dln_nifty','dln_india_ip','dln_usd_inr',     'dln_msci_world','d_us3mtbills', 2;
 'India C2+T (CPI)', 'dln_nifty','dln_india_ip','dln_india_cpi',   'dln_msci_world','d_us3mtbills', 2;
 'India C3+T (IR)',  'dln_nifty','dln_india_ip','d_indiacallrate', 'dln_msci_world','d_us3mtbills', 5;
 'India C4+T (M3)',  'dln_nifty','dln_india_ip','dln_india_m3',    'dln_msci_world','d_us3mtbills', 2;
 'Korea D1+T (FX)',  'dln_kospi','dln_korea_ip','dln_usd_kwr',     'dln_msci_world','d_us3mtbills', 2;
 'Korea D2+T (CPI)', 'dln_kospi','dln_korea_ip','dln_korea_cpi',   'dln_msci_world','d_us3mtbills', 2;
 'Korea D3+T (IR)',  'dln_kospi','dln_korea_ip','d_koreair',       'dln_msci_world','d_us3mtbills', 2;
 'Korea D4+T (M3)',  'dln_kospi','dln_korea_ip','dln_korea_m3',    'dln_msci_world','d_us3mtbills', 3;
};

%% ================================================================
%  BLOCK 1 — RESIDUAL DIAGNOSTICS (normality, ARCH, serial corr.)
% ================================================================
fprintf('================================================================\n');
fprintf('BLOCK 1: RESIDUAL DIAGNOSTICS (per equation, all 8 specs)\n');
fprintf('JB = Jarque-Bera normality (H0: normal)\n');
fprintf('ARCH = Engle ARCH test, 4 lags (H0: no ARCH)\n');
fprintf('LB = Ljung-Box, 12 lags (H0: no autocorrelation)\n');
fprintf('================================================================\n');

for i = 1:size(specs,1)
    name = specs{i,1}; vars = specs(i,2:6); p = specs{i,7};
    Y = data{:, vars}; Y = Y(all(~isnan(Y),2), :);
    Mdl = varm(5, p); EstMdl = estimate(Mdl, Y);
    E = infer(EstMdl, Y);   % residuals (T-p x 5)

    fprintf('\n--- %s (lag %d) ---\n', name, p);
    fprintf('%-16s %8s %8s %8s %8s\n','Equation','JB_p','ARCH_p','LB_p','skew/kurt');
    for k = 1:5
        e = E(:,k);
        [~, jb_p] = jbtest(e);
        [~, arch_p] = archtest(e, 'Lags', 4);
        [~, lb_p] = lbqtest(e, 'Lags', 12);
        fprintf('%-16s %8.3f %8.3f %8.3f   %.2f/%.1f\n', ...
            vars{k}, jb_p, arch_p, lb_p, skewness(e), kurtosis(e));
    end
end
fprintf('\n[Interpretation: expect JB rejected (non-normal) & some ARCH — this is\n');
fprintf(' exactly why bootstrap inference (Block 2) is preferable to asymptotic SEs.]\n');

%% ================================================================
%  BLOCK 2 — BOOTSTRAP FEVD OF IP + INDIA vs KOREA GAP TEST
%  Residual (iid) bootstrap: resample residuals, rebuild data, re-estimate,
%  recompute generalized FEVD of IP at month 12, normalize. B replications.
% ================================================================
fprintf('\n\n================================================================\n');
fprintf('BLOCK 2: BOOTSTRAP FEVD OF INDUSTRIAL PRODUCTION (month 12)\n');
fprintf('Generalized FEVD, row-normalized. 90%% percentile CIs.\n');
fprintf('================================================================\n');

B = 500;                 % replications (bump to 1000 for final run if time allows)
H = 12;                  % horizon
globalShareCI = struct();

for i = 1:size(specs,1)
    name = specs{i,1}; vars = specs(i,2:6); p = specs{i,7};
    Y = data{:, vars}; Y = Y(all(~isnan(Y),2), :);
    T = size(Y,1);
    ipIdx = 2; msciIdx = 4; tbillIdx = 5;   % positions in [stock,ip,dom,msci,tbill]

    Mdl = varm(5, p); EstMdl = estimate(Mdl, Y);
    E = infer(EstMdl, Y);                 % residuals
    Efull = [NaN(p,5); E];                % pad to align
    constC = EstMdl.Constant; AR = EstMdl.AR;

    bootGlobal = zeros(B,1);   % global (MSCI+Tbill) share of IP
    bootStock  = zeros(B,1);   % domestic stock share of IP
    for b = 1:B
        % resample residuals with replacement
        idx = randi(size(E,1), size(E,1), 1);
        Eb = E(idx,:);
        % rebuild series recursively from initial p obs
        Yb = zeros(T,5); Yb(1:p,:) = Y(1:p,:);
        for t = p+1:T
            xt = constC;
            for L = 1:p
                xt = xt + AR{L}*Yb(t-L,:)';
            end
            Yb(t,:) = xt' + Eb(t-p,:);
        end
        try
            Mb = estimate(varm(5,p), Yb);
            D = fevd(Mb, 'NumObs', H, 'Method','generalized');
            raw = squeeze(D(H,:,ipIdx));
            raw = raw / sum(raw);          % normalize
            bootGlobal(b) = raw(msciIdx) + raw(tbillIdx);
            bootStock(b)  = raw(1);
        catch
            bootGlobal(b) = NaN; bootStock(b) = NaN;
        end
    end
    bootGlobal = bootGlobal(~isnan(bootGlobal));
    g_lo = prctile(bootGlobal,5); g_hi = prctile(bootGlobal,95); g_pt = median(bootGlobal);

    fprintf('\n--- %s ---\n', name);
    fprintf('Global (MSCI+Tbill) share of IP: %.1f%%  [90%% CI: %.1f%%, %.1f%%]\n', ...
        100*g_pt, 100*g_lo, 100*g_hi);
    globalShareCI.(sprintf('spec%d',i)) = bootGlobal;
end

% ---- India vs Korea: is the global share significantly larger in India? ----
fprintf('\n=== INDIA vs KOREA GLOBAL-SHARE GAP (bootstrap test) ===\n');
india = [globalShareCI.spec1; globalShareCI.spec2; globalShareCI.spec4];   % drop C3 (not supported)
korea = [globalShareCI.spec5; globalShareCI.spec6; globalShareCI.spec7; globalShareCI.spec8];
% pooled difference distribution via random pairing
nrep = min(numel(india), numel(korea));
diffs = india(randperm(numel(india),nrep)) - korea(randperm(numel(korea),nrep));
fprintf('Mean India global share: %.1f%%\n', 100*mean(india));
fprintf('Mean Korea global share: %.1f%%\n', 100*mean(korea));
fprintf('India - Korea difference: %.1f pp  [90%% CI: %.1f, %.1f]\n', ...
    100*median(diffs), 100*prctile(diffs,5), 100*prctile(diffs,95));
fprintf('P(India > Korea) = %.3f\n', mean(diffs>0));
fprintf('[If the 90%% CI excludes 0 / P(India>Korea) near 1, the gap is significant,\n');
fprintf(' not just numerically larger — this directly strengthens H2 reversal claim.]\n');
