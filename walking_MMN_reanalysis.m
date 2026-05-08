function walking_MMN_reanalysis(data_dir, output_dir, ref_method, recursive)
%WALKING_MMN_REANALYSIS  Re-analysis of walking MMN data
%
% USAGE:
%   walking_MMN_reanalysis()
%       -> Opens dialogs to select input folder and create output folder.
%   walking_MMN_reanalysis(data_dir, output_dir, ref_method, recursive)
%       -> Programmatic mode (no dialogs).
%
% INPUTS:
%   data_dir   : folder with cs_new .mat files (omit -> dialog)
%   output_dir : where to save results (omit -> dialog with create option)
%   ref_method : 'common_avg' | 'original' | 'mastoid' | 'cms'
%                (omit -> dialog with selection list)
%   recursive  : true to search subdirectories (default: false)
%
% DIFFERENCES FROM ORIGINAL walking_MMN.m:
%   1. Reference: uses common-average of 32 EEG channels (CAR) by default,
%      configurable. Original used (X(:,37)+X(:,38)+X(:,39))/3 which had
%      high variance in cols 37, 38 - injecting noise.
%   2. Baseline correction: -50 to 0 ms mean subtracted (original had none).
%   3. dev-std convention (negative MMN, standard convention).
%      Original used std-dev (positive MMN).
%   4. Multiple time windows reported (P50, N1, classical MMN, late MMN, P3a)
%   5. Both mean and peak amplitudes reported per window.
%
% NOTE: ICA artifact removal is NOT included here. For matching cleanliness,
% post-process in Python with MNE or use EEGLAB ICA.

use_dialogs = (nargin < 1);

% === Step 1: Select INPUT folder (dialog if not specified) ===
if nargin < 1 || isempty(data_dir)
    data_dir = uigetdir(pwd, ...
        'Select folder containing .mat files (input)');
    if isequal(data_dir, 0)
        fprintf('Cancelled.\n');
        return;
    end
end
if ~exist(data_dir, 'dir')
    error('Input folder does not exist: %s', data_dir);
end

% === Step 2: Specify OUTPUT folder (with creation option) ===
if nargin < 2 || isempty(output_dir)
    parent_dir = uigetdir(data_dir, ...
        'Select PARENT folder for output (the new output folder will be created here)');
    if isequal(parent_dir, 0)
        fprintf('Cancelled.\n');
        return;
    end
    default_name = sprintf('reanalysis_%s', datestr(now, 'yyyymmdd_HHMMSS'));
    answer = inputdlg('Enter new output folder name:', ...
        'Output folder', 1, {default_name});
    if isempty(answer)
        fprintf('Cancelled.\n');
        return;
    end
    output_dir = fullfile(parent_dir, answer{1});
end
if ~exist(output_dir, 'dir')
    [ok, msg] = mkdir(output_dir);
    if ~ok
        error('Failed to create output folder: %s\n%s', output_dir, msg);
    end
    fprintf('Created output folder: %s\n', output_dir);
end

% === Step 3: Select reference method (dialog if not specified) ===
if nargin < 3 || isempty(ref_method)
    if use_dialogs
        ref_choices = {'common_avg', 'original', 'mastoid', 'cms'};
        ref_descs = {'Common Average Reference (recommended)', ...
                     'Original script: (X37+X38+X39)/3', ...
                     'Linked mastoid: (X37+X38)/2', ...
                     'No re-reference (Biosemi CMS)'};
        [sel, ok] = listdlg('PromptString', 'Select reference method:', ...
                            'SelectionMode', 'single', ...
                            'ListString', ref_descs, ...
                            'InitialValue', 1, ...
                            'ListSize', [400 100]);
        if ~ok
            fprintf('Cancelled.\n');
            return;
        end
        ref_method = ref_choices{sel};
    else
        ref_method = 'common_avg';
    end
end

if nargin < 4 || isempty(recursive)
    recursive = false;
end

% === Step 4: Find all .mat files ===
if recursive
    files = dir(fullfile(data_dir, '**', '*.mat'));
else
    files = dir(fullfile(data_dir, '*.mat'));
end
files = files(~[files.isdir]);
N = length(files);

if N == 0
    error('No .mat files found in %s (recursive=%d)', data_dir, recursive);
end

fprintf('\n========================================\n');
fprintf('walking_MMN_reanalysis\n');
fprintf('========================================\n');
fprintf('Input folder : %s\n', data_dir);
fprintf('Output folder: %s\n', output_dir);
fprintf('Reference    : %s\n', ref_method);
fprintf('Recursive    : %d\n', recursive);
fprintf('Files found  : %d\n', N);
fprintf('----------------------------------------\n');
for k = 1:N
    fprintf('  [%d] %s\n', k, files(k).name);
end
fprintf('========================================\n\n');

% === Step 5: Setup parameters ===
Fs = 1024;
n_pre  = round(50 * Fs/1000);
n_post = round(450 * Fs/1000);
epoch_len = n_pre + n_post + 1;
times_ms = ((0:epoch_len-1) - n_pre) / Fs * 1000;

ch_names = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5', ...
            'P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8', ...
            'CP6','CP2','C4','T8','FC6','FC2','F4','F8','AF4','Fp2', ...
            'Fz','Cz'};

windows = struct( ...
    'P50',     [40, 80], ...
    'N1',      [80, 130], ...
    'MMN_classic', [130, 200], ...
    'MMN_late',    [150, 250], ...
    'P3a',     [180, 280] );

% === Step 6: Pre-allocate storage ===
ave_std_all = zeros(epoch_len, 32, N);
ave_dev_all = zeros(epoch_len, 32, N);
n_std_all   = zeros(N, 1);
n_dev_all   = zeros(N, 1);
fnames      = cell(N, 1);

win_names = fieldnames(windows);
n_win = length(win_names);
features = zeros(N, 32 * n_win * 2 + 4);

% === Step 7: Main loop ===
% To use parallel processing, change 'for' to 'parfor' below
% (requires Parallel Computing Toolbox)
t0 = tic;
for k = 1:N
    fpath = fullfile(files(k).folder, files(k).name);
    fnames{k} = files(k).name;
    fprintf('\n[%d/%d] Processing %s\n', k, N, files(k).name);

    try
        [ave_std, ave_dev, n_std, n_dev] = process_one_subject( ...
            fpath, ref_method, Fs, n_pre, n_post);
    catch err
        warning('Failed to process %s: %s', files(k).name, err.message);
        continue;
    end

    ave_std_all(:,:,k) = ave_std;
    ave_dev_all(:,:,k) = ave_dev;
    n_std_all(k) = n_std;
    n_dev_all(k) = n_dev;
    fprintf('  n_std=%d, n_dev=%d\n', n_std, n_dev);

    diff_wave = ave_dev - ave_std;
    feat_row = zeros(1, 32 * n_win * 2 + 4);
    col = 1;
    for w = 1:n_win
        wlim = windows.(win_names{w});
        wmask = times_ms >= wlim(1) & times_ms <= wlim(2);
        for ch = 1:32
            feat_row(col) = mean(diff_wave(wmask, ch));
            feat_row(col + 1) = peak_in_window(diff_wave(wmask, ch));
            col = col + 2;
        end
    end
    feat_row(end-3:end) = [n_std, n_dev, n_std + n_dev, k];
    features(k, :) = feat_row;
end
elapsed = toc(t0);
fprintf('\nTotal processing time: %.1f sec\n', elapsed);

% === Step 8: Save outputs ===
write_features_csv(fullfile(output_dir, 'per_subject_features.csv'), ...
    features, fnames, ch_names, win_names);

for k = 1:N
    if n_std_all(k) == 0, continue; end
    [~, base, ~] = fileparts(fnames{k});
    csvwrite(fullfile(output_dir, sprintf('per_subject_aveStd_%s.csv', base)), ...
        ave_std_all(:,:,k));
    csvwrite(fullfile(output_dir, sprintf('per_subject_aveDev_%s.csv', base)), ...
        ave_dev_all(:,:,k));
end

valid_idx = n_std_all > 0;
ga_std  = mean(ave_std_all(:,:,valid_idx), 3);
ga_dev  = mean(ave_dev_all(:,:,valid_idx), 3);
ga_diff = ga_dev - ga_std;

csvwrite(fullfile(output_dir, 'grand_average_std.csv'),  ga_std);
csvwrite(fullfile(output_dir, 'grand_average_dev.csv'),  ga_dev);
csvwrite(fullfile(output_dir, 'grand_average_diff.csv'), ga_diff);
csvwrite(fullfile(output_dir, 'time_axis_ms.csv'), times_ms');

write_info(fullfile(output_dir, 'info.txt'), ...
    data_dir, output_dir, ref_method, N, Fs, n_pre, n_post, ...
    n_std_all, n_dev_all, fnames, win_names, windows, elapsed);

plot_grand_average(ga_std, ga_dev, ga_diff, times_ms, ...
    windows, output_dir, ref_method);

fprintf('\nDone. Results saved in:\n  %s\n', output_dir);
if use_dialogs
    msgbox(sprintf('Analysis complete.\n\nResults saved to:\n%s', output_dir), ...
        'walking_MMN_reanalysis');
end
end


function [ave_std, ave_dev, n_std_used, n_dev_used] = ...
    process_one_subject(fpath, ref_method, Fs, n_pre, n_post)

epoch_len = n_pre + n_post + 1;
loaded = load(fpath);
X = loaded.cs_new;
n_samp = size(X, 1);

switch ref_method
    case 'original'
        ave_ref = (X(:,37) + X(:,38) + X(:,39)) / 3;
    case 'mastoid'
        ave_ref = (X(:,37) + X(:,38)) / 2;
    case 'common_avg'
        ave_ref = mean(X(:, 1:32), 2);
    case 'cms'
        ave_ref = zeros(n_samp, 1);
    otherwise
        error('Unknown ref_method: %s', ref_method);
end

Y = zeros(n_samp, 36);
for i = 1:36
    Y(:, i) = X(:, i) - ave_ref;
end

[bn, an] = butter(4, [49 51] / (Fs/2), 'stop');
[bp, ap] = butter(4, [1 20] / (Fs/2), 'bandpass');
for i = 1:36
    Y(:, i) = filtfilt(bn, an, Y(:, i));
    Y(:, i) = filtfilt(bp, ap, Y(:, i));
end

EOG = Y(:, 35) - Y(:, 36);

[A, idx] = findpeaks(X(:, 41), 'MinPeakDistance', 100);
trig_all = [idx, A];
if trig_all(1, 1) <= 51
    trig_all = trig_all(2:end-1, :);
end

n_trig = size(trig_all, 1);
valid = true(n_trig, 1);
for t = 1:n_trig
    s = trig_all(t, 1) - n_pre;
    e = trig_all(t, 1) + n_post;
    if s < 1 || e > n_samp
        valid(t) = false; continue;
    end
    epoch_eeg = Y(s:e, 1:32);
    if any(abs(epoch_eeg(:)) > 100)
        valid(t) = false; continue;
    end
    epoch_eog = EOG(s:e);
    if any(abs(epoch_eog) > 100)
        valid(t) = false;
    end
end

trig_valid = trig_all(valid, :);
dev_mask = trig_valid(:,2) <= -2.106079e+05 & trig_valid(:,2) >= -2.1060792e+05;
std_mask = trig_valid(:,2) <= -2.1060793e+05;
dev_samples = trig_valid(dev_mask, 1);
std_samples = trig_valid(std_mask, 1);

[ave_dev, n_dev_used] = average_epochs(Y(:, 1:32), dev_samples, n_pre, n_post);
[ave_std, n_std_used] = average_epochs(Y(:, 1:32), std_samples, n_pre, n_post);

bl = 1:n_pre;
ave_dev = ave_dev - mean(ave_dev(bl, :), 1);
ave_std = ave_std - mean(ave_std(bl, :), 1);
end


function [ave, n_used] = average_epochs(Y, samples, n_pre, n_post)
epoch_len = n_pre + n_post + 1;
n_samp = size(Y, 1);
add = zeros(epoch_len, 32);
n_used = 0;
for j = 1:length(samples)
    s = samples(j) - n_pre;
    e = samples(j) + n_post;
    if s >= 1 && e <= n_samp
        add = add + Y(s:e, :);
        n_used = n_used + 1;
    end
end
if n_used > 0
    ave = add / n_used;
else
    ave = zeros(epoch_len, 32);
end
end


function val = peak_in_window(seg)
[pks_pos, ~] = findpeaks(seg);
[pks_neg, ~] = findpeaks(-seg);
if isempty(pks_pos) && isempty(pks_neg)
    val = 0; return;
end
if isempty(pks_neg)
    val = max(pks_pos);
elseif isempty(pks_pos)
    val = -max(pks_neg);
else
    pos_max = max(pks_pos);
    neg_max = max(pks_neg);
    if neg_max >= pos_max
        val = -neg_max;
    else
        val = pos_max;
    end
end
end


function write_features_csv(fpath, features, fnames, ch_names, win_names)
fid = fopen(fpath, 'w');
header = {'filename'};
for w = 1:length(win_names)
    for ch = 1:length(ch_names)
        header{end+1} = sprintf('%s_%s_mean', win_names{w}, ch_names{ch});
        header{end+1} = sprintf('%s_%s_peak', win_names{w}, ch_names{ch});
    end
end
header = [header, {'n_std','n_dev','n_total','subj_idx'}];
fprintf(fid, '%s\n', strjoin(header, ','));
for k = 1:size(features, 1)
    row = features(k, :);
    fprintf(fid, '%s', fnames{k});
    fprintf(fid, ',%g', row);
    fprintf(fid, '\n');
end
fclose(fid);
end


function write_info(fpath, data_dir, output_dir, ref_method, N, Fs, ...
    n_pre, n_post, n_std_all, n_dev_all, fnames, win_names, windows, elapsed)
fid = fopen(fpath, 'w');
fprintf(fid, 'walking_MMN_reanalysis run info\n');
fprintf(fid, '================================\n');
fprintf(fid, 'date: %s\n', datestr(now));
fprintf(fid, 'data_dir: %s\n', data_dir);
fprintf(fid, 'output_dir: %s\n', output_dir);
fprintf(fid, 'ref_method: %s\n', ref_method);
fprintf(fid, 'N files: %d\n', N);
fprintf(fid, 'Fs: %d Hz\n', Fs);
fprintf(fid, 'epoch: %d ms pre, %d ms post, total %d samples\n', ...
    round(n_pre*1000/Fs), round(n_post*1000/Fs), n_pre + n_post + 1);
fprintf(fid, 'baseline: -50 to 0 ms\n');
fprintf(fid, 'filters: 50 Hz notch (BW order 4), 1-20 Hz bandpass (BW order 4)\n');
fprintf(fid, 'rejection: |EEG ch1-32| > 100 uV OR |filtered EOG| > 100 uV\n');
fprintf(fid, 'elapsed: %.1f sec\n\n', elapsed);
fprintf(fid, 'time windows:\n');
for w = 1:length(win_names)
    wlim = windows.(win_names{w});
    fprintf(fid, '  %s: %d-%d ms\n', win_names{w}, wlim(1), wlim(2));
end
fprintf(fid, '\nper-file diagnostics:\n');
fprintf(fid, '%-50s %8s %8s\n', 'filename', 'n_std', 'n_dev');
for k = 1:N
    fprintf(fid, '%-50s %8d %8d\n', fnames{k}, n_std_all(k), n_dev_all(k));
end
fclose(fid);
end


function plot_grand_average(ga_std, ga_dev, ga_diff, times_ms, ...
    windows, output_dir, ref_method)
fig = figure('Position', [100 100 1200 500], 'Visible', 'off');
labels = {'Fz', 'Cz'};
chs    = [31, 32];
for c = 1:2
    subplot(1, 2, c); hold on;
    plot(times_ms, ga_std(:, chs(c)), 'k-', 'LineWidth', 1.5);
    plot(times_ms, ga_dev(:, chs(c)), 'r-', 'LineWidth', 1.5);
    plot(times_ms, ga_diff(:, chs(c)), 'b--', 'LineWidth', 2);
    yl = ylim;
    fn = fieldnames(windows);
    cmap = lines(length(fn));
    for w = 1:length(fn)
        wlim = windows.(fn{w});
        patch([wlim(1) wlim(2) wlim(2) wlim(1)], ...
              [yl(1) yl(1) yl(2) yl(2)], cmap(w,:), ...
              'EdgeColor', 'none', 'FaceAlpha', 0.10);
    end
    plot(times_ms, ga_std(:, chs(c)), 'k-', 'LineWidth', 1.5);
    plot(times_ms, ga_dev(:, chs(c)), 'r-', 'LineWidth', 1.5);
    plot(times_ms, ga_diff(:, chs(c)), 'b--', 'LineWidth', 2);
    set(gca, 'YDir', 'reverse');
    xlabel('Time (ms)'); ylabel('Amplitude (\muV)');
    legend({'Standard','Deviant','Diff (dev-std)'}, 'Location', 'best');
    title(sprintf('%s (ref: %s)', labels{c}, ref_method));
    xlim([-50, 400]); grid on;
end
saveas(fig, fullfile(output_dir, 'grand_average_FzCz.png'));
close(fig);
end