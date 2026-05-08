function rmanova_GG_n22(csv_path, output_dir)

if nargin < 1 || isempty(csv_path)
    csv_path = 'per_subject_n22_full.csv';
end
if nargin < 2 || isempty(output_dir)
    output_dir = pwd;
end

T = readtable(csv_path);
n = height(T);
fprintf('Loaded %d subjects from %s\n\n', n, csv_path);

measures = {'topo_sim', 'Fz_MMN', 'F4_MMN', 'lateralization'};
n_measures = length(measures);

WithinDesign = table({'walk';'water';'clay'}, ...
    'VariableNames', {'cond'});
WithinDesign.cond = categorical(WithinDesign.cond, ...
    {'walk','water','clay'}, 'Ordinal', false);

fprintf('================================================================\n');
fprintf('Analysis 2: 1-way RM-ANOVA (walk/water/clay), n=%d\n', n);
fprintf('Mauchly sphericity + Greenhouse-Geisser correction\n');
fprintf('================================================================\n');

results = cell(n_measures, 11);
parfor i = 1:n_measures
    m = measures{i};
    sub = T(:, {[ 'walk_' m], ['water_' m], ['clay_' m]});
    sub.Properties.VariableNames = {'walk','water','clay'};

    rm = fitrm(sub, 'walk-clay ~ 1', 'WithinDesign', WithinDesign);
    aov = ranova(rm);
    sphtbl = mauchly(rm);

    F = aov{'(Intercept):cond', 'F'};
    df1 = aov{'(Intercept):cond', 'DF'};
    df2 = aov{'Error(cond)', 'DF'};
    p_uncorr = aov{'(Intercept):cond', 'pValue'};
    p_gg = aov{'(Intercept):cond', 'pValueGG'};
    SS_eff = aov{'(Intercept):cond', 'SumSq'};
    SS_err = aov{'Error(cond)', 'SumSq'};
    eta2_p = SS_eff / (SS_eff + SS_err);

    W_sph = sphtbl.W;
    p_sph = sphtbl.pValue;

    eps_struct = epsilon(rm);
    eps_gg = eps_struct.GreenhouseGeisser;

    sph_ok = p_sph >= 0.05;

    walk_m = mean(sub.walk);
    water_m = mean(sub.water);
    clay_m = mean(sub.clay);

    results(i, :) = {m, walk_m, water_m, clay_m, ...
        W_sph, p_sph, sph_ok, eps_gg, F, p_uncorr, p_gg};
end

for i = 1:n_measures
    r = results(i, :);
    fprintf('\n--- %s ---\n', r{1});
    fprintf('  means: walk=%+.3f  water=%+.3f  clay=%+.3f\n', r{2}, r{3}, r{4});
    fprintf('  Mauchly: W=%.4f, p=%.4f  -> sphericity %s\n', ...
        r{5}, r{6}, ternary(r{7}, 'OK', 'VIOLATED'));
    fprintf('  GG epsilon = %.4f\n', r{8});
    fprintf('  Uncorrected:  F(2,42) = %.3f, p = %.4f\n', r{9}, r{10});
    fprintf('  GG-corrected: F(%.2f,%.2f) = %.3f, p = %.4f\n', ...
        2*r{8}, 42*r{8}, r{9}, r{11});
    if ~r{7}
        fprintf('  -> REPORT GG-CORRECTED p = %.4f\n', r{11});
    end
end

fprintf('\n================================================================\n');
fprintf('Bonferroni post-hoc (paired t, 3 comparisons, alpha=0.0167)\n');
fprintf('================================================================\n');
fprintf('%-16s %-13s %7s %4s %8s %8s %7s %4s\n', ...
    'measure', 'pair', 't', 'df', 'p_raw', 'p_Bonf', 'd_z', 'sig');

pairs = {{'walk','water'}, {'walk','clay'}, {'water','clay'}};
n_pairs = n_measures * length(pairs);
posthoc = cell(n_pairs, 8);
parfor k = 1:n_pairs
    [im, ip] = ind2sub_local(k, n_measures, length(pairs));
    m = measures{im};
    c1 = pairs{ip}{1}; c2 = pairs{ip}{2};
    x1 = T.([c1 '_' m]);
    x2 = T.([c2 '_' m]);
    diff = x1 - x2;
    [~, p_raw, ~, st] = ttest(x1, x2);
    p_bonf = min(p_raw * 3, 1);
    d = mean(diff) / std(diff);
    flag = '';
    if p_bonf < 0.01, flag = '**';
    elseif p_bonf < 0.05, flag = '*';
    end
    posthoc(k, :) = {m, [c1 '-' c2], st.tstat, st.df, p_raw, p_bonf, d, flag};
end

for i = 1:n_measures
    for j = 1:length(pairs)
        k = (i-1)*length(pairs) + j;
        r = posthoc(k, :);
        fprintf('%-16s %-13s %7.3f %4d %8.4f %8.4f %7.3f %4s\n', ...
            r{1}, r{2}, r{3}, r{4}, r{5}, r{6}, r{7}, r{8});
    end
end

A2 = cell2table(results, 'VariableNames', ...
    {'measure','walk_mean','water_mean','clay_mean', ...
     'mauchly_W','mauchly_p','sphericity_OK','GG_epsilon', ...
     'F','p_uncorrected','p_GG'});
PH = cell2table(posthoc, 'VariableNames', ...
    {'measure','pair','t','df','p_raw','p_Bonf','cohens_dz','sig'});

writetable(A2, fullfile(output_dir, 'results_analysis2_with_GG_matlab.csv'));
writetable(PH, fullfile(output_dir, 'results_posthoc_bonferroni_matlab.csv'));
fprintf('\nResults saved to %s\n', output_dir);

end


function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end

function [im, ij] = ind2sub_local(k, ~, n_inner)
im = ceil(k / n_inner);
ij = mod(k-1, n_inner) + 1;
end