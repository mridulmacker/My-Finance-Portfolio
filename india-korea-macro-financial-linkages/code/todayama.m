%% ================================================================
%  TABLE E3 — TODA-YAMAMOTO causality, exact p-values for headline arrows
%  Runs TY (VAR at p+1, test first p lags) and prints the specific
%  directional tests cited in Section 5.5.
% ================================================================
clear; clc;
load('master_dataset_transformed.mat');

specs = {
    'India FX',  'ln_nifty','ln_india_ip','ln_usd_inr',     'ln_msci_world','us3m_tbills', 2;
    'India CPI', 'ln_nifty','ln_india_ip','ln_india_cpi',   'ln_msci_world','us3m_tbills', 2;
    'India IR',  'ln_nifty','ln_india_ip','india_call_rate','ln_msci_world','us3m_tbills', 5;
    'India M3',  'ln_nifty','ln_india_ip','ln_india_m3',    'ln_msci_world','us3m_tbills', 2;
    'Korea FX',  'ln_kospi','ln_korea_ip','ln_usd_kwr',     'ln_msci_world','us3m_tbills', 2;
    'Korea CPI', 'ln_kospi','ln_korea_ip','ln_korea_cpi',   'ln_msci_world','us3m_tbills', 2;
    'Korea IR',  'ln_kospi','ln_korea_ip','korea_ir',       'ln_msci_world','us3m_tbills', 2;
    'Korea M3',  'ln_kospi','ln_korea_ip','ln_korea_m3',    'ln_msci_world','us3m_tbills', 3;
    };
% positions: 1=stock, 2=IP, 3=domestic, 4=MSCI, 5=Tbill

function pval = ty_test(Y, p, causeIdx, respIdx)
dmax = 1; k = p + dmax; n = size(Y,2);
Ylhs = Y(k+1:end,:); nobs = size(Ylhs,1);
X = ones(nobs,1);
for L=1:k, X=[X, Y(k+1-L:end-L,:)]; end
B = X\Ylhs; resid = Ylhs - X*B;
Sig = (resid'*resid)/(nobs-size(X,2));
idx = [];
for L=1:p, idx(end+1) = 1+(L-1)*n+causeIdx; end
b = B(idx,respIdx);
XtXinv = inv(X'*X);
V = Sig(respIdx,respIdx)*XtXinv(idx,idx);
W = b'*(V\b);
pval = 1 - chi2cdf(W, p);
end

fprintf('SPEC | Stock->IP | IP->Stock | MSCI->Stock | Tbill->IP | Tbill->dom(FXonly)\n');
fprintf('----------------------------------------------------------------------------\n');
for i=1:size(specs,1)
    name=specs{i,1}; vars=specs(i,2:6); p=specs{i,7};
    Y = data{:, vars}; Y = Y(all(~isnan(Y),2), :);
    p1 = ty_test(Y,p,1,2);   % stock -> IP
    p2 = ty_test(Y,p,2,1);   % IP -> stock
    p3 = ty_test(Y,p,4,1);   % MSCI -> stock
    p4 = ty_test(Y,p,5,2);   % Tbill -> IP
    p5 = ty_test(Y,p,5,3);   % Tbill -> domestic (won/rupee etc; meaningful in FX spec)
    fprintf('%-9s | %.4f | %.4f | %.4f | %.4f | %.4f\n', name, p1,p2,p3,p4,p5);
end
fprintf('\n(Stock=Nifty/KOSPI, IP=industrial production, MSCI=MSCI World, Tbill=US T-bill,\n');
fprintf(' dom=rotating domestic var. Tbill->dom is only economically meaningful in the FX spec.)\n');