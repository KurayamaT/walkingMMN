function sit_vs_walk_diagnostic()

sit_dir = uigetdir('', 'Select SIT (sitting) data folder containing *_sitting.mat files');
if isequal(sit_dir, 0); return; end

walk_dir = uigetdir(sit_dir, 'Select WALK (free) data folder containing *_free.mat files');
if isequal(walk_dir, 0); return; end

[parent_dir, ~] = fileparts(sit_dir);
out_input = inputdlg('Output folder name:', 'Output', 1, {'sit_vs_walk_diagnostic_out'});
if isempty(out_input); return; end
output_dir = fullfile(parent_dir, out_input{1});
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

fprintf('Sit dir   : %s\n', sit_dir);
fprintf('Walk dir  : %s\n', walk_dir);
fprintf('Output dir: %s\n\n', output_dir);

fs = 1024;
ch_names = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5', ...
            'P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8', ...
            'CP6','CP2','C4','T8','FC6','FC2','F4','F8','AF4','Fp2', ...
            'Fz','Cz'};
n_ch = 32;

[b_bp, a_bp] = butter(4, [1 20] / (fs/2), 'bandpass');

sit_files  = dir(fullfile(sit_dir,  '*_sitting.mat'));
walk_files = dir(fullfile(walk_dir, '*_free.mat'));

if isempty(sit_files)
    errordlg(sprintf('No *_sitting.mat files found in %s', sit_dir));
    return;
end
if isempty(walk_files)
    errordlg(sprintf('No *_free.mat files found in %s', walk_dir));
    return;
end

sit_subjects  = arrayfun(@(f) erase(f.name, '_sitting.mat'), sit_files,  'UniformOutput', false);
walk_subjects = arrayfun(@(f) erase(f.name, '_free.mat'),    walk_files, 'UniformOutput', false);
common = intersect(sit_subjects, walk_subjects);
n_sub = length(common);

if n_sub == 0
    errordlg('No common subjects found between Sit and Walk folders.');
    return;
end
fprintf('Found %d common subjects across Sit and Walk.\n\n', n_sub);

freqs_target = 0:0.5:50;
psd_sit  = zeros(n_sub, length(freqs_target), n_ch);
psd_walk = zeros(n_sub, length(freqs_target), n_ch);
trial_std_sit  = zeros(n_sub, 1);
trial_std_walk = zeros(n_sub, 1);

parfor si = 1:n_sub
    subj = common{si};
    fprintf('[%d/%d] %s\n', si, n_sub, subj);

    d_sit = load(fullfile(sit_dir, [subj '_sitting.mat']));
    X_sit = d_sit.cs_new(:, 1:n_ch);
    X_sit = X_sit - mean(X_sit, 2);
    X_sit = filtfilt(b_bp, a_bp, X_sit);
    [pxx_sit, ~] = pwelch(X_sit, fs*4, fs*2, freqs_target, fs);
    psd_sit(si, :, :) = pxx_sit;
    trial_std_sit(si) = mean(std(X_sit));

    d_walk = load(fullfile(walk_dir, [subj '_free.mat']));
    X_walk = d_walk.cs_new(:, 1:n_ch);
    X_walk = X_walk - mean(X_walk, 2);
    X_walk = filtfilt(b_bp, a_bp, X_walk);
    [pxx_walk, ~] = pwelch(X_walk, fs*4, fs*2, freqs_target, fs);
    psd_walk(si, :, :) = pxx_walk;
    trial_std_walk(si) = mean(std(X_walk));
end

mean_psd_sit  = squeeze(mean(psd_sit,  1));
mean_psd_walk = squeeze(mean(psd_walk, 1));

figure('Position', [100 100 1500 950], 'Name', 'Sit vs Walk diagnostic');

subplot(2, 3, 1);
target_chs = {'Fz','Cz','F4','Fp1','Fp2','Oz'};
for ci = 1:length(target_chs)
    idx = find(strcmp(ch_names, target_chs{ci}));
    semilogy(freqs_target, mean_psd_sit(:, idx), '-', 'LineWidth', 1.2,  ...
        'Color', [0 0.4 0.8], 'DisplayName', sprintf('Sit %s', target_chs{ci}));
    hold on;
end
for ci = 1:length(target_chs)
    idx = find(strcmp(ch_names, target_chs{ci}));
    semilogy(freqs_target, mean_psd_walk(:, idx), '--', 'LineWidth', 1.2, ...
        'Color', [0.85 0.3 0.1], 'DisplayName', sprintf('Walk %s', target_chs{ci}));
end
xlabel('Frequency (Hz)'); ylabel('PSD (\muV^2/Hz)');
title('PSD by channel: Sit (solid) vs Walk (dashed)');
legend('Location', 'eastoutside', 'FontSize', 7); grid on;
xlim([0 50]);

subplot(2, 3, 2);
psd_diff_log = 10*log10(mean_psd_walk ./ mean_psd_sit);
imagesc(freqs_target, 1:n_ch, psd_diff_log');
colormap(gca, redblue_cmap());
clim([-6 6]);
set(gca, 'YTick', 1:n_ch, 'YTickLabel', ch_names, 'FontSize', 7);
xlabel('Frequency (Hz)');
title('Walk - Sit (dB)');
cb = colorbar; ylabel(cb, 'dB');

subplot(2, 3, 3);
fz_idx = find(strcmp(ch_names, 'Fz'));
cz_idx = find(strcmp(ch_names, 'Cz'));
mean_diff_dB_fz = 10*log10(mean(psd_walk(:,:,fz_idx),1) ./ mean(psd_sit(:,:,fz_idx),1));
mean_diff_dB_cz = 10*log10(mean(psd_walk(:,:,cz_idx),1) ./ mean(psd_sit(:,:,cz_idx),1));
plot(freqs_target, mean_diff_dB_fz, 'b', 'LineWidth', 2); hold on;
plot(freqs_target, mean_diff_dB_cz, 'r', 'LineWidth', 2);
yline(0, 'k--');
xlabel('Frequency (Hz)'); ylabel('Walk-Sit (dB)');
legend('Fz', 'Cz', 'Location', 'best');
title('Spectral perturbation: Walk relative to Sit');
grid on; xlim([0 50]);

subplot(2, 3, 4);
plot([1, 2], [trial_std_sit, trial_std_walk]', 'Color', [0.7 0.7 0.7]); hold on;
plot([1, 2], [mean(trial_std_sit), mean(trial_std_walk)], 'k-', 'LineWidth', 3);
errorbar([1 2], [mean(trial_std_sit), mean(trial_std_walk)], ...
    [std(trial_std_sit)/sqrt(n_sub), std(trial_std_walk)/sqrt(n_sub)], 'k', 'LineWidth', 2);
[~, p_std] = ttest(trial_std_sit, trial_std_walk);
xticks([1 2]); xticklabels({'Sit','Walk'});
ylabel('Mean across-channel std (\muV)');
title(sprintf('Continuous std: paired t, p=%.4f', p_std));
grid on; xlim([0.5 2.5]);

subplot(2, 3, 5:6);
psd_band_diff = zeros(n_ch, 5);
band_labels = {'delta(1-3)', 'theta(4-7)', 'alpha(8-12)', 'beta(13-25)', 'low-gamma(25-40)'};
band_ranges = {[1 3], [4 7], [8 12], [13 25], [25 40]};
for bi = 1:5
    fmask = freqs_target >= band_ranges{bi}(1) & freqs_target <= band_ranges{bi}(2);
    band_sit  = squeeze(mean(mean_psd_sit(fmask,:), 1));
    band_walk = squeeze(mean(mean_psd_walk(fmask,:), 1));
    psd_band_diff(:, bi) = 10*log10(band_walk ./ band_sit);
end
bar(psd_band_diff);
set(gca, 'XTick', 1:n_ch, 'XTickLabel', ch_names, 'XTickLabelRotation', 90, 'FontSize', 7);
ylabel('Walk - Sit (dB)');
yline(0, 'k-');
legend(band_labels, 'Location', 'eastoutside');
title('Band-wise spectral perturbation per channel');
grid on;

print(gcf, fullfile(output_dir, 'sit_vs_walk_diagnostic.png'), '-dpng', '-r150');

T = table(common, trial_std_sit, trial_std_walk, ...
    'VariableNames', {'subject','sit_continuous_std','walk_continuous_std'});
writetable(T, fullfile(output_dir, 'sit_vs_walk_continuous_std.csv'));

writematrix(psd_band_diff, fullfile(output_dir, 'psd_band_diff_dB_per_channel.csv'));
writematrix(freqs_target',     fullfile(output_dir, 'psd_freqs.csv'));
writematrix(mean_psd_sit,      fullfile(output_dir, 'psd_sit_mean.csv'));
writematrix(mean_psd_walk,     fullfile(output_dir, 'psd_walk_mean.csv'));

fprintf('\nDone. Saved to %s\n', output_dir);

msgbox(sprintf(['Sit vs Walk diagnostic complete.\n\n' ...
    'n = %d common subjects\n' ...
    'Output: %s'], n_sub, output_dir), 'Done');

end


function cmap = redblue_cmap()
n = 64;
r = [linspace(0, 1, n/2), linspace(1, 1, n/2)]';
g = [linspace(0, 1, n/2), linspace(1, 0, n/2)]';
b = [linspace(1, 1, n/2), linspace(1, 0, n/2)]';
cmap = [r g b];
end