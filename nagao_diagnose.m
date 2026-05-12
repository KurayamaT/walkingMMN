clear; close all; clc;

[fname, fpath] = uigetfile({'nagao_*.mat;nagao*.mat', 'nagao files'; '*.mat', 'All MAT files'}, ...
    'nagao の .mat ファイルを選択（例: nagao_sitting.mat）');
if isequal(fname, 0); return; end

load(fullfile(fpath, fname));
fprintf('Loaded: %s\n', fname);
fprintf('cs_new size: [%d %d]\n\n', size(cs_new));

n_col = size(cs_new, 2);

fprintf('==================== 全列の概要 ====================\n');
fprintf('%4s  %14s %14s %14s %14s %10s\n', ...
    'Col', 'min', 'max', 'mean', 'std', 'unique_n');
for i = 1:n_col
    col = cs_new(:, i);
    n_uniq = length(unique(round(col, 2)));
    fprintf('%4d  %14.4f %14.4f %14.4f %14.4f %10d\n', ...
        i, min(col), max(col), mean(col), std(col), min(n_uniq, 99999));
end

fprintf('\n==================== トリガー列の候補 ====================\n');
fprintf('（離散的な値を持ち、特定値が多数回出現する列、または大振幅列）\n\n');

trigger_candidates = [];
for i = 33:n_col
    col = cs_new(:, i);
    col_round = round(col, 2);
    [vals, ~, ic] = unique(col_round);
    counts = accumarray(ic, 1);
    n_distinct = length(vals);
    [sorted_counts, sort_i] = sort(counts, 'descend');
    top5_pct = sum(sorted_counts(1:min(5, end))) / sum(counts) * 100;
    has_large_abs = any(abs(col) > 100000);
    
    if (n_distinct < 100 && top5_pct > 90) || has_large_abs
        trigger_candidates(end+1) = i;
        fprintf('Col %d: 離散値=%d, top5値=%.1f%%, 最大絶対値=%.1f', ...
            i, n_distinct, top5_pct, max(abs(col)));
        if has_large_abs
            fprintf(' [大振幅あり]');
        end
        fprintf('\n');
        for j = 1:min(5, length(vals))
            v = vals(sort_i(j));
            c = sorted_counts(j);
            fprintf('     値 %14.4f : %d 回 (%.2f%%)\n', ...
                v, c, c/sum(counts)*100);
        end
        fprintf('\n');
    end
end

if isempty(trigger_candidates)
    fprintf('離散値ベースで候補が見つかりませんでした。\n');
    fprintf('代わりに、列ごとの peak 数 を調べます...\n\n');
    for i = 33:n_col
        col = cs_new(:, i);
        [pk_amp, pk_idx] = findpeaks(col, 'MinPeakDistance', 100);
        if length(pk_idx) > 100 && length(pk_idx) < 5000
            fprintf('Col %d: %d peaks detected (range %.4f to %.4f)\n', ...
                i, length(pk_idx), min(col), max(col));
        end
        [npk_amp, ~] = findpeaks(-col, 'MinPeakDistance', 100);
        if length(npk_amp) > 100 && length(npk_amp) < 5000
            fprintf('Col %d (NEG): %d neg peaks detected\n', i, length(npk_amp));
        end
    end
end

fprintf('\n==================== Col 41 (標準形式) 詳細 ====================\n');
if n_col >= 41
    col41 = cs_new(:, 41);
    fprintf('範囲: %.4f to %.4f\n', min(col41), max(col41));
    fprintf('平均: %.4f, 標準偏差: %.4f\n', mean(col41), std(col41));
    [pk_amp_41, ~] = findpeaks(col41, 'MinPeakDistance', 100);
    fprintf('正方向 peak 数 (MinPeakDistance=100): %d\n', length(pk_amp_41));
    if ~isempty(pk_amp_41)
        edges = linspace(min(pk_amp_41), max(pk_amp_41), 11);
        hcounts = histcounts(pk_amp_41, edges);
        fprintf('Peak 振幅 10-bin histogram:\n');
        for k = 1:length(hcounts)
            fprintf('  [%10.4f, %10.4f): %d\n', edges(k), edges(k+1), hcounts(k));
        end
    end
end

fprintf('\n==================== Col 57 (imokawa 形式) 詳細 ====================\n');
if n_col >= 57
    col57 = cs_new(:, 57);
    fprintf('範囲: %.4f to %.4f\n', min(col57), max(col57));
    fprintf('平均: %.4f, 標準偏差: %.4f\n', mean(col57), std(col57));
    [pk_amp_57, ~] = findpeaks(col57, 'MinPeakDistance', 450);
    fprintf('正方向 peak 数 (MinPeakDistance=450): %d\n', length(pk_amp_57));
    if ~isempty(pk_amp_57)
        edges = linspace(min(pk_amp_57), max(pk_amp_57), 11);
        hcounts = histcounts(pk_amp_57, edges);
        fprintf('Peak 振幅 10-bin histogram:\n');
        for k = 1:length(hcounts)
            fprintf('  [%10.4f, %10.4f): %d\n', edges(k), edges(k+1), hcounts(k));
        end
    end
else
    fprintf('Col 57 は存在しません (n_col=%d)\n', n_col);
end

fs = 1024;
n_plot = min(30 * fs, size(cs_new, 1));

figure('Position', [50 50 1500 900], 'Name', sprintf('nagao_diagnose: %s', fname));

n_panels = 2 + length(trigger_candidates(1:min(3, end)));
if n_panels < 3; n_panels = 3; end

subplot(n_panels, 1, 1);
if n_col >= 41
    plot((1:n_plot)/fs, cs_new(1:n_plot, 41));
    title('Col 41 (standard-format trigger)');
else
    text(0.5, 0.5, 'Col 41 not present', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Time (s)'); ylabel('Value'); grid on;

subplot(n_panels, 1, 2);
if n_col >= 57
    plot((1:n_plot)/fs, cs_new(1:n_plot, 57));
    title('Col 57 (imokawa-format trigger)');
else
    text(0.5, 0.5, 'Col 57 not present', 'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Time (s)'); ylabel('Value'); grid on;

for k = 1:min(length(trigger_candidates), n_panels - 2)
    subplot(n_panels, 1, 2 + k);
    cc = trigger_candidates(k);
    plot((1:n_plot)/fs, cs_new(1:n_plot, cc));
    title(sprintf('Col %d (candidate %d)', cc, k));
    xlabel('Time (s)'); ylabel('Value'); grid on;
end

[~, base, ~] = fileparts(fname);
saveas(gcf, fullfile(fpath, sprintf('%s_diagnose.png', base)));

fprintf('\n');
fprintf('====================================================\n');
fprintf('SUMMARY (これを Claude に貼り付けてください)\n');
fprintf('====================================================\n');
fprintf('File: %s\n', fname);
fprintf('cs_new size: [%d %d]\n', size(cs_new));
if n_col >= 41
    fprintf('Col 41 range: %.4f to %.4f (std=%.4f)\n', min(col41), max(col41), std(col41));
end
if n_col >= 57
    fprintf('Col 57 range: %.4f to %.4f (std=%.4f)\n', min(col57), max(col57), std(col57));
end
fprintf('Trigger candidates: %s\n', mat2str(trigger_candidates));
if ~isempty(trigger_candidates)
    fprintf('Most likely trigger column: %d\n', trigger_candidates(1));
end
fprintf('Figure saved: %s_diagnose.png\n', base);
fprintf('====================================================\n');
