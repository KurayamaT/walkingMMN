
function walking_MMN_imokawa()
% walking_MMN_imokawa.m  (CORRECTED based on user's original 1778221030855_imokawa_walking_MMN.m)
%
% imokawa subject specific re-analysis.
%
% IMPORTANT — verified from user's original script:
%   - Trigger column is Col 57 (NOT Col 41 as in standard format)
%   - cs_new is 624640 x 57
%   - Std trigger: -161535.583 <= value <= -161535.56
%   - Dev trigger: value > -161535.56
%   - MinPeakDistance = 450 samples (~440 ms)
%
% Output: same format as walking_MMN_singletrial.m for direct integration
%   into the n=21 group (becoming n=22).
%
% Run once per condition (sitting / free / walter / nenndo).

% ==================== GUI ====================

input_dir = uigetdir('', 'Select directory containing imokawa .mat files');
if isequal(input_dir, 0); return; end

[parent_dir, default_name] = fileparts(input_dir);
out_input = inputdlg('Output folder name:', 'Output', 1, ...
    {[default_name '_imokawa']});
if isempty(out_input); return; end
output_dir = fullfile(parent_dir, out_input{1});
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

ref_options = {'Common Average Reference (CAR, recommended)', ...
               'Original mastoid-like (X37+X38+X39)/3'};
[ref_idx, ok] = listdlg('PromptString', 'Reference method:', ...
                        'SelectionMode', 'single', ...
                        'ListString', ref_options, ...
                        'ListSize', [400, 100]);
if ~ok; return; end

cond_input = inputdlg('Condition suffix (sitting / free / walter / nenndo):', ...
    'Condition', 1, {'sitting'});
if isempty(cond_input); return; end
cond_label = cond_input{1};

% ==================== Constants ====================

fs = 1024;
n_pre = 51;
n_post = 460;
epoch_len = n_pre + n_post + 1;
times_ms = ((0:epoch_len-1) - n_pre) / fs * 1000;
reject_uv = 100;

ch_names = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5', ...
            'P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8', ...
            'CP6','CP2','C4','T8','FC6','FC2','F4','F8','AF4','Fp2', ...
            'Fz','Cz'};

[b_bp, a_bp] = butter(4, [1 20] / (fs/2), 'bandpass');
[b_th, a_th] = butter(4, [4 8] / (fs/2), 'bandpass');

mmn_mask = times_ms >= 130 & times_ms <= 200;
bl_mask  = times_ms <  0;

left_chs  = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5'};
right_chs = {'Fp2','AF4','F8','F4','FC2','FC6','T8','C4','CP2','CP6'};
[~, l_idx] = ismember(left_chs,  ch_names);
[~, r_idx] = ismember(right_chs, ch_names);

% ==================== File discovery ====================

target_pattern = ['imokawa_' cond_label '.mat'];
target_path = fullfile(input_dir, target_pattern);

if ~exist(target_path, 'file')
    files = dir(fullfile(input_dir, 'imokawa*.mat'));
    target_file = '';
    for fi = 1:length(files)
        if strcmpi(files(fi).name, target_pattern)
            target_file = files(fi).name;
            target_path = fullfile(input_dir, target_file);
            break;
        end
    end
    if isempty(target_file)
        errordlg(sprintf('File not found: %s', target_pattern));
        return;
    end
end

fprintf('Processing: %s (condition=%s)\n\n', target_pattern, cond_label);

% ==================== Load and process ====================

data = load(target_path);
if ~isfield(data, 'cs_new')
    errordlg('cs_new not found in file');
    return;
end
X = data.cs_new;
fprintf('Loaded cs_new: [%d x %d]\n', size(X));

if size(X, 2) < 57
    errordlg(sprintf('Need at least 57 columns. Got %d.', size(X, 2)));
    return;
end

% --- Re-reference ---
switch ref_idx
    case 1, ref_signal = mean(X(:, 1:32), 2);
    case 2, ref_signal = mean(X(:, [37 38 39]), 2);
end
eeg = X(:, 1:32) - ref_signal;
fprintf('Re-referenced (%s)\n', ref_options{ref_idx});

% --- Filter ---
eeg_filt = filtfilt(b_bp, a_bp, eeg);
fprintf('Bandpass-filtered 1-20 Hz\n');

% --- Trigger detection (imokawa-specific: Col 57) ---
trig = X(:, 57);
fprintf('\nTrigger detection (imokawa format, Col 57):\n');
fprintf('  Col 57 range: %.4f to %.4f\n', min(trig), max(trig));

[pk_amp, pk_idx] = findpeaks(trig, 'MinPeakDistance', 450);
fprintf('  Peaks detected (MinPeakDistance=450): %d\n', length(pk_amp));

% Original script's classification:
% std: -161535.583 <= val <= -161535.56
% dev: val > -161535.56
std_mask = pk_amp >= -161535.583 & pk_amp <= -161535.56;
dev_mask = pk_amp > -161535.56;

std_idx = pk_idx(std_mask);
dev_idx = pk_idx(dev_mask);
n_std_raw = length(std_idx);
n_dev_raw = length(dev_idx);

fprintf('  Initial counts: std=%d, dev=%d\n', n_std_raw, n_dev_raw);
if n_std_raw + n_dev_raw > 0
    fprintf('  Std/Dev ratio: %.1f%% / %.1f%%\n', ...
        n_std_raw/(n_std_raw+n_dev_raw)*100, ...
        n_dev_raw/(n_std_raw+n_dev_raw)*100);
end

if ~isempty(std_idx) && std_idx(1) <= n_pre
    std_idx = std_idx(2:end);
end
if ~isempty(dev_idx) && dev_idx(1) <= n_pre
    dev_idx = dev_idx(2:end);
end

% --- Epoch + reject ---
[epochs_std, k_std] = make_epochs(eeg_filt, std_idx, n_pre, n_post, reject_uv);
[epochs_dev, k_dev] = make_epochs(eeg_filt, dev_idx, n_pre, n_post, reject_uv);
n_std = sum(k_std);
n_dev = sum(k_dev);

fprintf('\n  After ±%g µV rejection:\n', reject_uv);
fprintf('    std: %d → %d (kept %.1f%%)\n', n_std_raw, n_std, n_std/n_std_raw*100);
fprintf('    dev: %d → %d (kept %.1f%%)\n', n_dev_raw, n_dev, n_dev/n_dev_raw*100);

if n_std < 10 || n_dev < 10
    errordlg(sprintf('Insufficient trials. std=%d, dev=%d', n_std, n_dev));
    return;
end

% --- Baseline correct ---
epochs_std = baseline_correct(epochs_std, bl_mask);
epochs_dev = baseline_correct(epochs_dev, bl_mask);

% --- Output ---
out_subdir = fullfile(output_dir, sprintf('imokawa_%s', cond_label));
if ~exist(out_subdir, 'dir'); mkdir(out_subdir); end

ave_std = mean(epochs_std, 3);
ave_dev = mean(epochs_dev, 3);
diff_wave = ave_dev - ave_std;

writematrix(ave_std,   fullfile(out_subdir, 'aveStd.csv'));
writematrix(ave_dev,   fullfile(out_subdir, 'aveDev.csv'));
writematrix(diff_wave, fullfile(out_subdir, 'diffWave.csv'));

% ==================== TIER 1 indicators ====================

mmn_topo = mean(diff_wave(mmn_mask, :), 1);
lateralization = mean(mmn_topo(r_idx)) - mean(mmn_topo(l_idx));
gfp = std(mmn_topo);

fz_idx = find(strcmp(ch_names, 'Fz'));
cz_idx = find(strcmp(ch_names, 'Cz'));
f4_idx = find(strcmp(ch_names, 'F4'));

search_mask = times_ms >= 100 & times_ms <= 250;
search_t = times_ms(search_mask);
[~, peak_i] = min(diff_wave(search_mask, fz_idx));
peak_lat_fz = search_t(peak_i);

fz_mmn = mean(diff_wave(mmn_mask, fz_idx));
cz_mmn = mean(diff_wave(mmn_mask, cz_idx));
f4_mmn = mean(diff_wave(mmn_mask, f4_idx));

% ==================== TIER 2 single-trial indicators ====================

n_t = epoch_len;
n_ch = 32;

itc_dev = zeros(n_t, n_ch);
itc_std = zeros(n_t, n_ch);
theta_power_dev = zeros(n_t, n_ch);
theta_power_std = zeros(n_t, n_ch);
sd_dev = zeros(n_t, n_ch);
sd_std = zeros(n_t, n_ch);

epochs_dev_loc = epochs_dev;
epochs_std_loc = epochs_std;

parfor ci = 1:n_ch
    ch_dev = squeeze(epochs_dev_loc(:, ci, :));
    ch_std = squeeze(epochs_std_loc(:, ci, :));
    
    ch_dev_theta = filtfilt(b_th, a_th, ch_dev);
    ch_std_theta = filtfilt(b_th, a_th, ch_std);
    
    h_dev = hilbert(ch_dev_theta);
    h_std = hilbert(ch_std_theta);
    
    theta_power_dev(:, ci) = mean(abs(h_dev).^2, 2);
    theta_power_std(:, ci) = mean(abs(h_std).^2, 2);
    
    phase_dev = angle(h_dev);
    phase_std = angle(h_std);
    itc_dev(:, ci) = abs(mean(exp(1i * phase_dev), 2));
    itc_std(:, ci) = abs(mean(exp(1i * phase_std), 2));
    
    sd_dev(:, ci) = std(ch_dev, 0, 2);
    sd_std(:, ci) = std(ch_std, 0, 2);
end

writematrix(itc_dev,         fullfile(out_subdir, 'itc_dev.csv'));
writematrix(itc_std,         fullfile(out_subdir, 'itc_std.csv'));
writematrix(theta_power_dev, fullfile(out_subdir, 'theta_power_dev.csv'));
writematrix(theta_power_std, fullfile(out_subdir, 'theta_power_std.csv'));
writematrix(sd_dev,          fullfile(out_subdir, 'singletrial_sd_dev.csv'));
writematrix(sd_std,          fullfile(out_subdir, 'singletrial_sd_std.csv'));

bh_t = zeros(n_ch, 1);
bh_p = ones(n_ch, 1);
for ci = 1:n_ch
    dev_amp = squeeze(mean(epochs_dev(mmn_mask, ci, :), 1));
    std_amp = squeeze(mean(epochs_std(mmn_mask, ci, :), 1));
    [~, p_val, ~, st] = ttest2(dev_amp, std_amp);
    bh_t(ci) = st.tstat;
    bh_p(ci) = p_val;
end
bh_table = table(ch_names', bh_t, bh_p, ...
    'VariableNames', {'Channel', 'tstat', 'p'});
writetable(bh_table, fullfile(out_subdir, 'bishop_hardiman.csv'));

if n_dev >= 20
    half = floor(n_dev / 2);
    half1 = mean(epochs_dev(:, :, 1:half), 3) - ave_std;
    half2 = mean(epochs_dev(:, :, half+1:2*half), 3) - ave_std;
    topo1 = mean(half1(mmn_mask, :), 1);
    topo2 = mean(half2(mmn_mask, :), 1);
    sh_r = corrcoef(topo1, topo2);
    split_half_r = sh_r(1, 2);
else
    split_half_r = NaN;
end

ind_table = table({'imokawa'}, {cond_label}, n_std, n_dev, ...
    lateralization, gfp, peak_lat_fz, fz_mmn, cz_mmn, f4_mmn, split_half_r, ...
    'VariableNames', {'subject', 'condition', 'n_std', 'n_dev', ...
    'lateralization_RminusL_uV', 'GFP_uV', 'peak_lat_Fz_ms', ...
    'Fz_MMN_130_200_uV', 'Cz_MMN_130_200_uV', 'F4_MMN_130_200_uV', ...
    'split_half_topo_r'});
writetable(ind_table, fullfile(out_subdir, 'topo_indicators.csv'));

group_csv = fullfile(output_dir, 'group_summary.csv');
if exist(group_csv, 'file')
    try
        existing = readtable(group_csv);
        new_vars = ind_table.Properties.VariableNames;
        for vi = 1:length(new_vars)
            v = new_vars{vi};
            if any(strcmp(existing.Properties.VariableNames, v))
                if iscell(existing.(v)) && isnumeric(ind_table.(v))
                    existing.(v) = cellfun(@(x) str2double(string(x)), existing.(v));
                end
            end
        end
        common_cols = intersect(existing.Properties.VariableNames, new_vars);
        if length(common_cols) == width(existing)
            ind_subset = ind_table(:, common_cols);
            ind_subset.Properties.VariableNames = existing.Properties.VariableNames;
            combined = [existing; ind_subset];
            writetable(combined, group_csv);
        else
            writetable(ind_table, group_csv);
        end
    catch ME
        warning('Could not append to %s (%s). Writing fresh file.', ...
            group_csv, ME.message);
        writetable(ind_table, group_csv);
    end
else
    writetable(ind_table, group_csv);
end

writematrix(times_ms', fullfile(output_dir, 'time_axis_ms.csv'));

% ==================== Visualization ====================

figure('Position', [100 100 1400 800], 'Name', sprintf('imokawa_%s', cond_label));

ch_lbls = {'Fz','Cz','F4'};
for col_i = 1:3
    ci_plot = find(strcmp(ch_names, ch_lbls{col_i}));
    subplot(3, 3, col_i);
    plot(times_ms, ave_std(:, ci_plot), 'b', 'LineWidth', 1.5); hold on;
    plot(times_ms, ave_dev(:, ci_plot), 'r', 'LineWidth', 1.5);
    plot(times_ms, diff_wave(:, ci_plot), 'k', 'LineWidth', 2);
    yl_curr = get(gca, 'YLim');
    fill([130 200 200 130], [yl_curr(1) yl_curr(1) yl_curr(2) yl_curr(2)], ...
        [0 1 0], 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    xlim([-50 350]); set(gca, 'YDir', 'reverse');
    title(sprintf('imokawa %s (%s)', cond_label, ch_lbls{col_i}));
    xlabel('Time (ms)'); ylabel('µV');
    if col_i == 1
        legend('Std','Dev','Diff','MMN window', 'Location', 'best');
    end
    grid on;
end

subplot(3, 3, 4);
histogram(pk_amp, 50);
hold on;
yl = ylim;
plot([-161535.56 -161535.56], yl, 'r--', 'LineWidth', 2);
title('Trigger amplitude (Col 57)');
xlabel('Value'); ylabel('Count');
legend('Peaks', 'std/dev threshold', 'Location', 'best');
grid on;

subplot(3, 3, 5);
plot(times_ms, diff_wave); xlim([-50 350]); set(gca, 'YDir', 'reverse');
title('Difference wave (32 channels)');
xlabel('Time (ms)'); ylabel('µV'); grid on;

subplot(3, 3, 6);
bar(1:32, mmn_topo);
set(gca, 'XTickLabel', ch_names, 'XTickLabelRotation', 90, 'XTick', 1:32);
title('MMN topography (130-200 ms)');
ylabel('µV'); grid on;

subplot(3, 3, 7:9);
sample_n = min(30 * fs, length(trig));
plot((1:sample_n)/fs, trig(1:sample_n), 'b'); hold on;
std_in = std_idx(std_idx < sample_n);
dev_in = dev_idx(dev_idx < sample_n);
plot(std_in/fs, trig(std_in), 'go', 'MarkerSize', 8, 'LineWidth', 1.5);
plot(dev_in/fs, trig(dev_in), 'r^', 'MarkerSize', 10, 'LineWidth', 1.5);
title(sprintf('Trigger detection (first 30 sec): %d std, %d dev visible', ...
    length(std_in), length(dev_in)));
xlabel('Time (s)'); ylabel('Col 57 value');
legend('Col 57', 'Std', 'Dev');
grid on;

saveas(gcf, fullfile(out_subdir, 'imokawa_diagnostic.png'));

fprintf('\n');
fprintf('================================================\n');
fprintf('SUMMARY: imokawa %s\n', cond_label);
fprintf('================================================\n');
fprintf('  n_std kept: %d (initially %d)\n', n_std, n_std_raw);
fprintf('  n_dev kept: %d (initially %d)\n', n_dev, n_dev_raw);
fprintf('  Fz MMN (130-200 ms):  %+.2f µV\n', fz_mmn);
fprintf('  Cz MMN (130-200 ms):  %+.2f µV\n', cz_mmn);
fprintf('  F4 MMN (130-200 ms):  %+.2f µV\n', f4_mmn);
fprintf('  Peak latency (Fz):    %.0f ms\n', peak_lat_fz);
fprintf('  Lateralization (R-L): %+.3f µV\n', lateralization);
fprintf('  GFP:                  %.3f µV\n', gfp);
fprintf('  Split-half topo r:    %.3f\n', split_half_r);
fprintf('  Output:               %s\n', out_subdir);
fprintf('================================================\n');
fprintf('Reference typical values from n=21:\n');
fprintf('  Fz MMN sit:    -0.68 µV (group mean)\n');
fprintf('  Topo split-r:  >0.5 = good\n');
fprintf('================================================\n');

msgbox(sprintf(['imokawa_%s done.\n\n' ...
    'n_std=%d, n_dev=%d\n' ...
    'Fz MMN: %+.2f µV (peak %.0f ms)\n\n' ...
    'Output: %s'], cond_label, n_std, n_dev, fz_mmn, peak_lat_fz, out_subdir), ...
    'Done');

end


% ==================== Helper functions ====================

function [epochs, keep] = make_epochs(eeg, idx_list, n_pre, n_post, reject_uv)
n_total = length(idx_list);
epoch_len = n_pre + n_post + 1;
n_ch = size(eeg, 2);
epochs_all = zeros(epoch_len, n_ch, n_total);
keep = false(n_total, 1);
for k = 1:n_total
    s = idx_list(k) - n_pre;
    e = idx_list(k) + n_post;
    if s < 1 || e > size(eeg, 1); continue; end
    ep = eeg(s:e, :);
    if any(abs(ep(:)) > reject_uv); continue; end
    epochs_all(:, :, k) = ep;
    keep(k) = true;
end
epochs = epochs_all(:, :, keep);
end

function epochs_corr = baseline_correct(epochs, bl_mask)
n_trials = size(epochs, 3);
epochs_corr = epochs;
for k = 1:n_trials
    bl_mean = mean(epochs(bl_mask, :, k), 1);
    epochs_corr(:, :, k) = epochs(:, :, k) - bl_mean;
end
end
