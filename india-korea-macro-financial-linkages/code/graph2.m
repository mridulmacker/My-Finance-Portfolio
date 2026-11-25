%% ================================================================
%  FIGURE 2 — Generalised IRF panels showing the channel difference
%  India IP response to a T-bill shock; Korea KOSPI response to MSCI shock.
% ================================================================
clear; clc;
load('master_dataset_transformed.mat');

H = 12;

% --- India FX spec: IP response to T-bill shock ---
viI = {'ln_nifty','ln_india_ip','ln_usd_inr','ln_msci_world','us3m_tbills'};
YI = data{:, viI}; YI = YI(all(~isnan(YI),2),:);
MI = estimate(varm(5,2), YI);
irfI = irf(MI, 'NumObs', H, 'Method','generalized');
% irf dims: (horizon, response, shock). Response=IP(2), Shock=Tbill(5)
respI = squeeze(irfI(:,2,5));

% --- Korea FX spec: KOSPI response to MSCI shock ---
viK = {'ln_kospi','ln_korea_ip','ln_usd_kwr','ln_msci_world','us3m_tbills'};
YK = data{:, viK}; YK = YK(all(~isnan(YK),2),:);
MK = estimate(varm(5,2), YK);
irfK = irf(MK, 'NumObs', H, 'Method','generalized');
% Response=KOSPI(1), Shock=MSCI(4)
respK = squeeze(irfK(:,1,4));

figure('Position',[100 100 900 380]);
subplot(1,2,1);
plot(0:H-1, respI, '-o','LineWidth',1.2); grid on;
title('(a) India: IP response to a US T-bill shock');
xlabel('Months'); ylabel('Response'); yline(0,'k-');

subplot(1,2,2);
plot(0:H-1, respK, '-o','LineWidth',1.2); grid on;
title('(b) Korea: KOSPI response to an MSCI World shock');
xlabel('Months'); ylabel('Response'); yline(0,'k-');

saveas(gcf, 'Figure2_IRF.png');
fprintf('Saved Figure2_IRF.png\n');