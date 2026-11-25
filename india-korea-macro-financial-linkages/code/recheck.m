%% ================================================================
%  TABLE 5.1 RE-RUN — CORRECTED: date-based COVID dummies and
%  pre-pandemic mask, computed fresh for EACH spec after that spec's
%  own missing-value filtering (fixes the Korea M3 misalignment and
%  the 1-row drift in the differenced specs).
% ================================================================
clear; clc;
load('master_dataset_transformed.mat');

% Build (or confirm) a full calendar date vector aligned to the RAW table
% BEFORE any filtering. Adjust the start date if your series doesn't
% begin January 2000 — check data.Date(1) or equivalent if you have one.
if ismember('Date', data.Properties.VariableNames)
    fullDates = data.Date;
elseif ismember('date', data.Properties.VariableNames)
    fullDates = data.date;
else
    warning('No date column found — constructing assumed calendar Jan-2000 start. VERIFY this matches your actual first observation before trusting results.');
    fullDates = (datetime(2000,1,1):calmonths(1):datetime(2000,1,1)+calmonths(height(data)-1))';
end

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

H = 12;
covidStart = datetime(2020,3,1); covidEnd = datetime(2020,5,31);

fprintf('=========================================================\n');
fprintf('PART A: FEVD of IP WITH COVID impulse dummies (date-based)\n');
fprintf('=========================================================\n');
for i = 1:size(specs,1)
    name=specs{i,1}; vars=specs(i,2:6); p=specs{i,7};
    Yraw = data{:, vars};
    validRows = all(~isnan(Yraw),2);
    Y = Yraw(validRows,:);
    datesThis = fullDates(validRows);          % <-- dates filtered THE SAME WAY as Y

    covidMask = (datesThis >= covidStart) & (datesThis <= covidEnd);
    fprintf('\n--- %s: %d COVID rows detected at dates: ', name, sum(covidMask));
    fprintf('%s ', datestr(datesThis(covidMask),'mmm-yyyy')); fprintf('---\n');
    if sum(covidMask) ~= 3
        warning('Expected 3 COVID months, found %d for %s — check date alignment before trusting this spec.', sum(covidMask), name);
    end

    D = double(covidMask);   % single dummy column (or split into 3 impulse dummies as before if preferred)
    T = size(Y,1);
    Dfull = zeros(T, sum(covidMask));
    idx = find(covidMask);
    for j=1:numel(idx), Dfull(idx(j),j)=1; end

    Mdl = varm(5,p);
    EstMdl = estimate(Mdl, Y, 'X', Dfull);
    Dc = fevd(EstMdl, 'NumObs', H, 'Method','generalized');
    raw = squeeze(Dc(H,:,2));
    normv = raw/sum(raw);
    fprintf('  Domestic stock: %.1f%%  MSCI: %.1f%%  T-bill: %.1f%%  Global total: %.1f%%\n', ...
        100*normv(1), 100*normv(4), 100*normv(5), 100*(normv(4)+normv(5)));
end

fprintf('\n\n=========================================================\n');
fprintf('PART B: FEVD of IP, PRE-PANDEMIC SUBSAMPLE (date-based, < Mar 2020)\n');
fprintf('=========================================================\n');
for i = 1:size(specs,1)
    name=specs{i,1}; vars=specs(i,2:6); p=specs{i,7};
    Yraw = data{:, vars};
    validRows = all(~isnan(Yraw),2);
    Y = Yraw(validRows,:);
    datesThis = fullDates(validRows);

    preMask = datesThis < covidStart;
    Ypre = Y(preMask,:);
    fprintf('\n--- %s: pre-pandemic sample = %d obs, ending %s ---\n', ...
        name, sum(preMask), datestr(datesThis(find(preMask,1,'last')),'mmm-yyyy'));

    Mdl = varm(5,p);
    EstMdl = estimate(Mdl, Ypre);
    Dc = fevd(EstMdl, 'NumObs', H, 'Method','generalized');
    raw = squeeze(Dc(H,:,2));
    normv = raw/sum(raw);
    fprintf('  Domestic stock: %.1f%%  MSCI: %.1f%%  T-bill: %.1f%%  Global total: %.1f%%\n', ...
        100*normv(1), 100*normv(4), 100*normv(5), 100*(normv(4)+normv(5)));
end