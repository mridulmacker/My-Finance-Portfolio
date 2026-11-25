%% ================================================================
%  FIGURE 4 — Rolling 60-month global-factor FEVD share, India vs Korea
%  Re-estimates the FX spec on rolling windows; plots global total share.
% ================================================================

clear; clc;
load('master_dataset_transformed.mat');
win = 60; H = 12;

specsC = {
    'India','dln_nifty','dln_india_ip','dln_usd_inr','dln_msci_world','d_us3mtbills',2;
    'Korea','dln_kospi','dln_korea_ip','dln_usd_kwr','dln_msci_world','d_us3mtbills',2;};

figure('Position',[100 100 850 450]); hold on;
colors = {'b','r'};
for s=1:2
    vars=specsC(s,2:6); p=specsC{s,7};
    Yall=data{:,vars}; valid=all(~isnan(Yall),2); Yall=Yall(valid,:);
    T=size(Yall,1);
    shares=nan(T,1);
    for t=win:T
        Yw=Yall(t-win+1:t,:);
        try
            M=estimate(varm(5,p),Yw);
            D=fevd(M,'NumObs',H,'Method','generalized');
            raw=squeeze(D(H,:,2)); raw=raw/sum(raw);
            shares(t)=100*(raw(4)+raw(5));   % MSCI + Tbill
        catch, end
    end
    plot(win:T, shares(win:T), colors{s}, 'LineWidth',1.3);
end
grid on; xlabel('Window end (month index)'); ylabel('Global-factor share of IP variance (%)');
legend('India','South Korea','Location','best');
title('Rolling 60-month global-factor share of industrial-production variance');
saveas(gcf,'Figure4_rolling.png');
fprintf('Saved Figure4_rolling.png\n');