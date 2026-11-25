load('master_dataset_transformed.mat');
x = data.dln_india_m3; x = x(~isnan(x));
for lag = [4 8 12]
    [h,p,stat] = kpsstest(x,'lags',lag);
    fprintf('KPSS M3 growth, %2d lags: stat=%.4f, p=%.4f -> %s\n', lag, stat, p, string(h==0));
end