%% ================================================================
%  E3 EXTRA COLUMNS — MSCI World -> IP, and domestic macro -> stock
%  (fixed: temp variable for the inverse before indexing)
% ================================================================
clear; clc;
load('master_dataset_transformed.mat');
specs = {
    'India FX','ln_nifty','ln_india_ip','ln_usd_inr','ln_msci_world','us3m_tbills',2;
    'India CPI','ln_nifty','ln_india_ip','ln_india_cpi','ln_msci_world','us3m_tbills',2;
    'India IR','ln_nifty','ln_india_ip','india_call_rate','ln_msci_world','us3m_tbills',5;
    'India M3','ln_nifty','ln_india_ip','ln_india_m3','ln_msci_world','us3m_tbills',2;
    'Korea FX','ln_kospi','ln_korea_ip','ln_usd_kwr','ln_msci_world','us3m_tbills',2;
    'Korea CPI','ln_kospi','ln_korea_ip','ln_korea_cpi','ln_msci_world','us3m_tbills',2;
    'Korea IR','ln_kospi','ln_korea_ip','korea_ir','ln_msci_world','us3m_tbills',2;
    'Korea M3','ln_kospi','ln_korea_ip','ln_korea_m3','ln_msci_world','us3m_tbills',3;};

function pval = ty(Y,p,cause,resp)
dmax=1; k=p+dmax; n=size(Y,2);
Ylhs=Y(k+1:end,:); nobs=size(Ylhs,1);
X=ones(nobs,1); for L=1:k, X=[X, Y(k+1-L:end-L,:)]; end
B=X\Ylhs; res=Ylhs-X*B; Sig=(res'*res)/(nobs-size(X,2));
idx=[]; for L=1:p, idx(end+1)=1+(L-1)*n+cause; end
b=B(idx,resp);
XtXinv = inv(X'*X);          % <-- temp variable fixes the syntax error
V = Sig(resp,resp)*XtXinv(idx,idx);
pval = 1 - chi2cdf(b'*(V\b), p);
end

fprintf('SPEC | MSCI->IP | domMacro->Stock\n');
fprintf('-----------------------------------\n');
for i=1:8
    vars=specs(i,2:6); p=specs{i,7};
    Y=data{:,vars}; Y=Y(all(~isnan(Y),2),:);
    p_msci_ip = ty(Y,p,4,2);   % MSCI World -> IP
    p_dom_stk = ty(Y,p,3,1);   % domestic macro -> stock
    fprintf('%-9s | %.4f | %.4f\n', specs{i,1}, p_msci_ip, p_dom_stk);
end