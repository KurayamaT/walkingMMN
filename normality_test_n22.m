function normality_test_n22(csv_path, output_dir)

if nargin < 1 || isempty(csv_path)
    csv_path = 'per_subject_n22_full.csv';
end
if nargin < 2 || isempty(output_dir)
    output_dir = pwd;
end

if ~exist(csv_path, 'file')
    error('CSV not found: %s', csv_path);
end
if ~exist(output_dir, 'dir')
    mkdir(output_dir);
end

T = readtable(csv_path);
n = height(T);
fprintf('Loaded %d subjects from %s\n\n', n, csv_path);

measures = {'topo_sim', 'Fz_MMN', 'F4_MMN', 'lateralization'};
n_measures = length(measures);

fprintf('================================================================\n');
fprintf('ANALYSIS 1: Sit - Walk paired differences (Shapiro-Wilk)\n');
fprintf('================================================================\n');
fprintf('%-18s %8s %8s %7s %7s %9s %9s %5s\n', ...
    'measure', 'mean', 'SD', 'skew', 'kurt', 'W', 'p', 'flag');

a1_results = cell(n_measures, 7);
parfor i = 1:n_measures
    m = measures{i};
    diff = T.(['sit_' m]) - T.(['walk_' m]);
    [W, p] = swtest_inline(diff);
    a1_results(i, :) = {m, mean(diff), std(diff), ...
        skewness(diff), kurtosis(diff)-3, W, p};
end
for i = 1:n_measures
    r = a1_results(i, :);
    fprintf('%-18s %8.3f %8.3f %7.2f %7.2f %9.4f %9.4f %5s\n', ...
        r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, sig_flag(r{7}));
end
fprintf('\n');

fprintf('================================================================\n');
fprintf('ANALYSIS 2: per cell (Walk, Water, Clay) Shapiro-Wilk\n');
fprintf('================================================================\n');
fprintf('%-18s %-7s %7s %7s %9s %9s %5s\n', ...
    'measure', 'cond', 'skew', 'kurt', 'W', 'p', 'flag');

conds_a2 = {'walk', 'water', 'clay'};
n_a2_cells = n_measures * 3;
a2_cell_results = cell(n_a2_cells, 6);
parfor k = 1:n_a2_cells
    [im, ic] = ind2sub_local(k, n_measures, 3);
    m = measures{im};
    c = conds_a2{ic};
    x = T.([c '_' m]);
    [W, p] = swtest_inline(x);
    a2_cell_results(k, :) = {m, c, skewness(x), kurtosis(x)-3, W, p};
end
for i = 1:n_measures
    for j = 1:3
        k = (i-1)*3 + j;
        r = a2_cell_results(k, :);
        fprintf('%-18s %-7s %7.2f %7.2f %9.4f %9.4f %5s\n', ...
            r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, sig_flag(r{6}));
    end
    fprintf('\n');
end

fprintf('================================================================\n');
fprintf('ANALYSIS 2: pairwise differences (post-hoc用)\n');
fprintf('================================================================\n');
fprintf('%-18s %-13s %7s %7s %9s %9s %5s\n', ...
    'measure', 'pair', 'skew', 'kurt', 'W', 'p', 'flag');

pairs = {{'walk','water'}, {'walk','clay'}, {'water','clay'}};
n_a2_pairs = n_measures * length(pairs);
a2_pair_results = cell(n_a2_pairs, 6);
parfor k = 1:n_a2_pairs
    [im, ip] = ind2sub_local(k, n_measures, length(pairs));
    m = measures{im};
    c1 = pairs{ip}{1}; c2 = pairs{ip}{2};
    diff = T.([c1 '_' m]) - T.([c2 '_' m]);
    [W, p] = swtest_inline(diff);
    a2_pair_results(k, :) = {m, [c1 '-' c2], ...
        skewness(diff), kurtosis(diff)-3, W, p};
end
for i = 1:n_measures
    for j = 1:length(pairs)
        k = (i-1)*length(pairs) + j;
        r = a2_pair_results(k, :);
        fprintf('%-18s %-13s %7.2f %7.2f %9.4f %9.4f %5s\n', ...
            r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, sig_flag(r{6}));
    end
    fprintf('\n');
end

A1 = cell2table(a1_results, 'VariableNames', ...
    {'measure','mean_diff','sd_diff','skew','kurt_excess','W','p'});
A2C = cell2table(a2_cell_results, 'VariableNames', ...
    {'measure','cond','skew','kurt_excess','W','p'});
A2P = cell2table(a2_pair_results, 'VariableNames', ...
    {'measure','pair','skew','kurt_excess','W','p'});

writetable(A1,  fullfile(output_dir, 'normality_analysis1.csv'));
writetable(A2C, fullfile(output_dir, 'normality_analysis2_cells.csv'));
writetable(A2P, fullfile(output_dir, 'normality_analysis2_pairs.csv'));

fprintf('================================================================\n');
fprintf('SUMMARY\n');
fprintf('================================================================\n');
n_v_a1  = sum(cell2mat(a1_results(:,7)) < 0.05);
n_v_a2c = sum(cell2mat(a2_cell_results(:,6)) < 0.05);
n_v_a2p = sum(cell2mat(a2_pair_results(:,6)) < 0.05);
fprintf('Analysis 1 violations  (p<.05): %d / %d\n', n_v_a1, n_measures);
fprintf('Analysis 2 cells       (p<.05): %d / %d\n', n_v_a2c, n_a2_cells);
fprintf('Analysis 2 pairs       (p<.05): %d / %d\n', n_v_a2p, n_a2_pairs);
fprintf('\nResults saved to %s\n', output_dir);

end


function [im, ij] = ind2sub_local(k, ~, n_inner)
im = ceil(k / n_inner);
ij = mod(k-1, n_inner) + 1;
end

function s = sig_flag(p)
if p < 0.001, s = '***';
elseif p < 0.01, s = '**';
elseif p < 0.05, s = '*';
else, s = '';
end
end

function [W, pValue] = swtest_inline(x)
x = x(~isnan(x));
x = sort(x(:));
n = length(x);
if n < 3
    W = NaN; pValue = NaN; return;
end

mtilde = norminv(((1:n)' - 3/8) ./ (n + 1/4));
m2 = sum(mtilde.^2);

u = 1/sqrt(n);
a = zeros(n, 1);
a(n) = -2.706056*u^5 + 4.434685*u^4 - 2.071190*u^3 ...
       - 0.147981*u^2 + 0.221157*u + mtilde(n)/sqrt(m2);
a(1) = -a(n);

if n > 5
    a(n-1) = -3.582633*u^5 + 5.682633*u^4 - 1.752460*u^3 ...
             - 0.293762*u^2 + 0.042981*u + mtilde(n-1)/sqrt(m2);
    a(2) = -a(n-1);
    eps_val = (m2 - 2*mtilde(n)^2 - 2*mtilde(n-1)^2) / ...
              (1 - 2*a(n)^2 - 2*a(n-1)^2);
    a(3:n-2) = mtilde(3:n-2) / sqrt(eps_val);
else
    eps_val = (m2 - 2*mtilde(n)^2) / (1 - 2*a(n)^2);
    a(2:n-1) = mtilde(2:n-1) / sqrt(eps_val);
end

mu_x = mean(x);
SS = sum((x - mu_x).^2);
if SS <= 0
    W = NaN; pValue = NaN; return;
end
W = (sum(a .* x))^2 / SS;

if n <= 11
    gamma_v = -2.273 + 0.459*n;
    mu = 0.5440 - 0.39978*n + 0.025054*n^2 - 0.0006714*n^3;
    sigma = exp(1.3822 - 0.77857*n + 0.062767*n^2 - 0.0020322*n^3);
    w = -log(gamma_v - log(1 - W));
    pValue = 1 - normcdf((w - mu) / sigma);
else
    lnN = log(n);
    mu = -1.5861 - 0.31082*lnN - 0.083751*lnN^2 + 0.0038915*lnN^3;
    sigma = exp(-0.4803 - 0.082676*lnN + 0.0030302*lnN^2);
    w = log(1 - W);
    pValue = 1 - normcdf((w - mu) / sigma);
end
end