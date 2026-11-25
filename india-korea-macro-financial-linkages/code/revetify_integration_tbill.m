%% ================================================================
%  RE-VERIFY INTEGRATION ORDER — T-bill and India call rate
%  Using both ADF (H0: unit root) and KPSS (H0: stationary) for
%  corroborating evidence, per Ben's feedback.
%  ================================================================
clear; clc;
load('master_dataset_transformed.mat');

varsToCheck = {'us3m_tbills', 'india_call_rate', 'korea_ir'};

for i = 1:numel(varsToCheck)
    v = varsToCheck{i};
    x = data.(v);
    x = x(~isnan(x));

    fprintf('\n=== %s ===\n', v);

    % ADF at multiple lag counts (H0: unit root / non-stationary)
    for lag = [4, 8, 12]
        [h_adf, p_adf, stat_adf] = adftest(x, 'lags', lag, 'model', 'ARD');
        verdict_adf = "Non-stationary";
        if h_adf == 1
            verdict_adf = "Stationary";
        end
        fprintf('  ADF, %2d lags: stat=%.4f, p=%.4f -> %s\n', lag, stat_adf, p_adf, verdict_adf);
    end

    % KPSS (H0: stationary -- OPPOSITE null, useful corroboration)
    for lag = [4, 8, 12]
        [h_kpss, p_kpss, stat_kpss] = kpsstest(x, 'lags', lag);
        verdict_kpss = "Stationary (fail to reject H0)";
        if h_kpss == 1
            verdict_kpss = "Non-stationary (reject H0)";
        end
        fprintf('  KPSS,%2d lags: stat=%.4f, p=%.4f -> %s\n', lag, stat_kpss, p_kpss, verdict_kpss);
    end
end

fprintf('\n\nINTERPRETATION GUIDE:\n');
fprintf('- ADF says stationary AND KPSS says stationary -> confidently I(0)\n');
fprintf('- ADF says non-stationary AND KPSS says non-stationary -> confidently I(1)\n');
fprintf('- Results disagree -> genuinely ambiguous, exactly Ben''s concern\n');