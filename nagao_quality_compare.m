function nagao_quality_compare()

base_dir = uigetdir('', 'Select base directory (parent of sitting_mat/free_mat/water_mat/nenndo_mat)');
if isequal(base_dir, 0); return; end

out_dir = fullfile(base_dir, 'nagao_quality_comparison');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

ref_options = {'Common Average Reference (CAR, recommended)', ...
               'Original mastoid-like (X37+X38+X39)/3'};
[ref_idx, ok] = listdlg('PromptString', 'Reference method:', ...
                        'SelectionMode', 'single', ...
                        'ListString', ref_options, ...
                        'ListSize', [400 100]);
if ~ok; return; end

conditions = {'sitting', 'free', 'water', 'nenndo'};
cond_dirs  = {'sitting_mat', 'free_mat', 'water_mat', 'nenndo_mat'};
subjects   = {'nagao', 'nagao2'};

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

n1_mask  = times_ms >= 80  & times_ms <= 130;
mmn_mask = times_ms >= 130 & times_ms <= 200;
bl_mask  = times_ms <  0;

fz = find(strcmp(ch_names, 'Fz'));
cz = find(strcmp(ch_names, 'Cz'));
f4 = find(strcmp(ch_names, 'F4'));

n_cond = length(conditions);
n_subj = length(subjects);

results = cell(n_cond * n_subj, 1);
diffwaves = cell(n_cond, n_subj);
ave_std_all = cell(n_cond, n_subj);
ave_dev_all = cell(n_cond, n_subj);
format_used = cell(n_cond, n_subj);

row = 0;
for ci = 1:n_cond
    cond = conditions{ci};
    cdir = fullfile(base_dir, cond_dirs{ci});
    if ~exist(cdir, 'dir')
        warning('Directory not found: %s', cdir);
        continue;
    end
    for si = 1:n_subj
        subj = subjects{si};
        fpath = fullfile(cdir, sprintf('%s_%s.mat', subj, cond));
        row = row + 1;
        
        if ~exist(fpath, 'file')
            fprintf('[%d/%d] MISSING: %s_%s.mat\n', row, n_cond*n_subj, subj, cond);
            results{row} = make_missing_row(subj, cond);
            diffwaves{ci, si} = NaN(epoch_len, 32);
            continue;
        end
        
        fprintf('[%d/%d] Processing %s_%s.mat ...\n', row, n_cond*n_subj, subj, cond);
        
        try
            tmp = load(fpath);
        catch ME
            fprintf('  ERROR: load failed: %s\n', ME.message);
            results{row} = make_error_row(subj, cond, 'load_failed');
            diffwaves{ci, si} = NaN(epoch_len, 32);
            continue;
        end
        
        if ~isfield(tmp, 'cs_new')
            fprintf('  ERROR: cs_new field not found\n');
            results{row} = make_error_row(subj, cond, 'no_cs_new');
            diffwaves{ci, si} = NaN(epoch_len, 32);
            continue;
        end
        
        X = tmp.cs_new;
        [n_samp, n_col] = size(X);
        
        switch ref_idx
            case 1, ref_signal = mean(X(:, 1:32), 2);
            case 2
                if n_col >= 39
                    ref_signal = mean(X(:, [37 38 39]), 2);
                else
                    ref_signal = mean(X(:, 1:32), 2);
                end
        end
        eeg = X(:, 1:32) - ref_signal;
        eeg_filt = filtfilt(b_bp, a_bp, eeg);
        
        [trig_col, trig_format, std_idx, dev_idx, pk_amp] = ...
            detect_format_and_extract_triggers(X, n_col);
        
        format_used{ci, si} = trig_format;
        
        if isempty(std_idx) || isempty(dev_idx)
            fprintf('  ERROR: no triggers detected (format=%s)\n', trig_format);
            results{row} = make_error_row(subj, cond, 'no_triggers');
            diffwaves{ci, si} = NaN(epoch_len, 32);
            continue;
        end
        
        n_std_raw = length(std_idx);
        n_dev_raw = length(dev_idx);
        
        if ~isempty(std_idx) && std_idx(1) <= n_pre
            std_idx = std_idx(2:end);
        end
        if ~isempty(dev_idx) && dev_idx(1) <= n_pre
            dev_idx = dev_idx(2:end);
        end
        
        [epochs_std, k_std] = make_epochs(eeg_filt, std_idx, n_pre, n_post, reject_uv);
        [epochs_dev, k_dev] = make_epochs(eeg_filt, dev_idx, n_pre, n_post, reject_uv);
        n_std = sum(k_std);
        n_dev = sum(k_dev);
        rej_rate_std = 1 - n_std / max(n_std_raw, 1);
        rej_rate_dev = 1 - n_dev / max(n_dev_raw, 1);
        
        if n_std < 10 || n_dev < 10
            fprintf('  WARN: insufficient trials (std=%d dev=%d)\n', n_std, n_dev);
            results{row} = make_insufficient_row(subj, cond, n_std, n_dev, ...
                n_std_raw, n_dev_raw, trig_col, trig_format);
            diffwaves{ci, si} = NaN(epoch_len, 32);
            continue;
        end
        
        epochs_std = baseline_correct(epochs_std, bl_mask);
        epochs_dev = baseline_correct(epochs_dev, bl_mask);
        
        ave_std = mean(epochs_std, 3);
        ave_dev = mean(epochs_dev, 3);
        diff_wave = ave_dev - ave_std;
        diffwaves{ci, si} = diff_wave;
        ave_std_all{ci, si} = ave_std;
        ave_dev_all{ci, si} = ave_dev;
        
        bl_sd_dev = median(std(epochs_dev(bl_mask, :, :), 0, 1), 'all');
        bl_sd_std = median(std(epochs_std(bl_mask, :, :), 0, 1), 'all');
        post_sd_dev = median(std(epochs_dev(~bl_mask, :, :), 0, 1), 'all');
        
        n1_fz = mean(diff_wave(n1_mask, fz));
        n1_cz = mean(diff_wave(n1_mask, cz));
        n1_f4 = mean(diff_wave(n1_mask, f4));
        
        mmn_fz = mean(diff_wave(mmn_mask, fz));
        mmn_cz = mean(diff_wave(mmn_mask, cz));
        mmn_f4 = mean(diff_wave(mmn_mask, f4));
        
        search_mask = times_ms >= 100 & times_ms <= 250;
        search_t = times_ms(search_mask);
        wave_fz = diff_wave(search_mask, fz);
        [peak_amp_fz, peak_i] = min(wave_fz);
        peak_lat_fz = search_t(peak_i);
        
        snr_fz = abs(peak_amp_fz) / max(bl_sd_dev, 1e-9);
        
        mmn_topo = mean(diff_wave(mmn_mask, :), 1);
        gfp = std(mmn_topo);
        
        half = floor(n_dev / 2);
        if half >= 5
            h1 = mean(epochs_dev(:, :, 1:half), 3) - ave_std;
            h2 = mean(epochs_dev(:, :, half+1:2*half), 3) - ave_std;
            t1 = mean(h1(mmn_mask, :), 1);
            t2 = mean(h2(mmn_mask, :), 1);
            r_mat = corrcoef(t1, t2);
            split_half_r = r_mat(1, 2);
        else
            split_half_r = NaN;
        end
        
        bh_p = ones(32, 1);
        bh_t = zeros(32, 1);
        epochs_dev_loc = epochs_dev;
        epochs_std_loc = epochs_std;
        parfor cc = 1:32
            dev_amp = squeeze(mean(epochs_dev_loc(mmn_mask, cc, :), 1));
            std_amp = squeeze(mean(epochs_std_loc(mmn_mask, cc, :), 1));
            [~, p_val, ~, st] = ttest2(dev_amp, std_amp);
            bh_t(cc) = st.tstat;
            bh_p(cc) = p_val;
        end
        n_sig_p05  = sum(bh_p < 0.05);
        n_sig_p01  = sum(bh_p < 0.01);
        n_sig_p001 = sum(bh_p < 0.001);
        
        results{row} = table({subj}, {cond}, {trig_format}, trig_col, ...
            n_std_raw, n_dev_raw, n_std, n_dev, rej_rate_std, rej_rate_dev, ...
            bl_sd_dev, bl_sd_std, post_sd_dev, ...
            n1_fz, n1_cz, n1_f4, ...
            mmn_fz, mmn_cz, mmn_f4, ...
            peak_amp_fz, peak_lat_fz, snr_fz, ...
            gfp, split_half_r, ...
            n_sig_p05, n_sig_p01, n_sig_p001, ...
            'VariableNames', {'subject','condition','format','trig_col', ...
            'n_std_raw','n_dev_raw','n_std_kept','n_dev_kept', ...
            'rej_rate_std','rej_rate_dev', ...
            'bl_sd_dev_uV','bl_sd_std_uV','post_sd_dev_uV', ...
            'N1_Fz_uV','N1_Cz_uV','N1_F4_uV', ...
            'MMN_Fz_uV','MMN_Cz_uV','MMN_F4_uV', ...
            'peak_amp_Fz_uV','peak_lat_Fz_ms','SNR_Fz', ...
            'GFP_uV','split_half_topo_r', ...
            'n_sig_ch_p05','n_sig_ch_p01','n_sig_ch_p001'});
        
        fprintf('  format=%s | n_dev=%d kept=%d | Fz MMN=%+.2f µV peak_lat=%.0fms SNR=%.2f | split_r=%.3f | BH(p<.05)=%d ch\n', ...
            trig_format, n_dev_raw, n_dev, mmn_fz, peak_lat_fz, snr_fz, split_half_r, n_sig_p05);
    end
end

valid_rows = ~cellfun(@isempty, results);
T = vertcat(results{valid_rows});
writetable(T, fullfile(out_dir, 'nagao_quality_comparison_full.csv'));
fprintf('\nWritten: nagao_quality_comparison_full.csv (%d rows)\n', height(T));

print_comparison_report(T, conditions, subjects, out_dir);

plot_diagnostic_panel(diffwaves, ave_std_all, ave_dev_all, ...
    conditions, subjects, times_ms, ch_names, fz, cz, f4, ...
    n1_mask, mmn_mask, out_dir);

fprintf('\nDone. Output: %s\n', out_dir);
msgbox(sprintf(['Quality comparison done.\n\nOutput: %s\n\n' ...
    'Inspect:\n  - nagao_quality_comparison_full.csv\n' ...
    '  - nagao_quality_comparison_report.txt\n' ...
    '  - nagao_quality_diagnostic.png'], out_dir), 'Done');

end


function [trig_col, trig_format, std_idx, dev_idx, pk_amp] = ...
    detect_format_and_extract_triggers(X, n_col)

trig_col = NaN;
trig_format = 'unknown';
std_idx = [];
dev_idx = [];
pk_amp = [];

if n_col >= 57
    trig_57 = X(:, 57);
    range_57 = max(trig_57) - min(trig_57);
    
    if any(trig_57 > -161536) && any(trig_57 < -161536 + 1) && range_57 > 0.01
        [pk_amp, pk_idx] = findpeaks(trig_57, 'MinPeakDistance', 450);
        std_mask = pk_amp >= -161535.583 & pk_amp <= -161535.56;
        dev_mask = pk_amp > -161535.56;
        std_idx = pk_idx(std_mask);
        dev_idx = pk_idx(dev_mask);
        trig_col = 57;
        trig_format = 'imokawa-like (Col 57)';
        return;
    end
end

if n_col >= 41
    trig_41 = X(:, 41);
    [pk_amp, pk_idx] = findpeaks(trig_41, 'MinPeakDistance', 100);
    std_mask_41 = pk_amp <= -210607.93;
    dev_mask_41 = pk_amp >= -210607.92 & pk_amp <= -210607.9;
    if sum(std_mask_41) > 100 && sum(dev_mask_41) > 10
        std_idx = pk_idx(std_mask_41);
        dev_idx = pk_idx(dev_mask_41);
        trig_col = 41;
        trig_format = 'standard (Col 41)';
        return;
    end
end

end


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


function row = make_missing_row(subj, cond)
row = table({subj}, {cond}, {'MISSING'}, NaN, ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
    'VariableNames', {'subject','condition','format','trig_col', ...
    'n_std_raw','n_dev_raw','n_std_kept','n_dev_kept', ...
    'rej_rate_std','rej_rate_dev', ...
    'bl_sd_dev_uV','bl_sd_std_uV','post_sd_dev_uV', ...
    'N1_Fz_uV','N1_Cz_uV','N1_F4_uV', ...
    'MMN_Fz_uV','MMN_Cz_uV','MMN_F4_uV', ...
    'peak_amp_Fz_uV','peak_lat_Fz_ms','SNR_Fz', ...
    'GFP_uV','split_half_topo_r', ...
    'n_sig_ch_p05','n_sig_ch_p01','n_sig_ch_p001'});
end


function row = make_error_row(subj, cond, reason)
row = make_missing_row(subj, cond);
row.format{1} = reason;
end


function row = make_insufficient_row(subj, cond, n_std, n_dev, n_std_raw, n_dev_raw, trig_col, trig_format)
row = make_missing_row(subj, cond);
row.format{1} = sprintf('insufficient (%s)', trig_format);
row.trig_col = trig_col;
row.n_std_raw = n_std_raw;
row.n_dev_raw = n_dev_raw;
row.n_std_kept = n_std;
row.n_dev_kept = n_dev;
end


function print_comparison_report(T, conditions, subjects, out_dir)

fid = fopen(fullfile(out_dir, 'nagao_quality_comparison_report.txt'), 'w');
print_to(fid, '====================================================================');
print_to(fid, '   NAGAO vs NAGAO2  QUALITY COMPARISON REPORT');
print_to(fid, '   Generated: %s', datestr(now));
print_to(fid, '====================================================================');
print_to(fid, '');
print_to(fid, 'Quality metrics (higher = better unless marked LOW=better):');
print_to(fid, '  - n_dev_kept:        retained deviant trials (HIGH=better)');
print_to(fid, '  - rej_rate_dev:      epoch rejection rate (LOW=better)');
print_to(fid, '  - bl_sd_dev:         baseline SD across deviant epochs (LOW=better)');
print_to(fid, '  - |N1_Fz|:           N1 amplitude at Fz (HIGH=better)');
print_to(fid, '  - |MMN_Fz|:          MMN-window amplitude at Fz (HIGH=better)');
print_to(fid, '  - SNR_Fz:            |peak amp| / baseline SD (HIGH=better)');
print_to(fid, '  - split_half_topo_r: within-subject reproducibility (HIGH=better)');
print_to(fid, '  - n_sig_ch_p05:      Bishop-Hardiman significant channel count (HIGH=better)');
print_to(fid, '');

wins = struct('nagao', 0, 'nagao2', 0, 'tie', 0);

for ci = 1:length(conditions)
    cond = conditions{ci};
    rows_c = strcmp(T.condition, cond);
    T_c = T(rows_c, :);
    
    print_to(fid, '');
    print_to(fid, '--------------------------------------------------------------------');
    print_to(fid, '  CONDITION: %s', upper(cond));
    print_to(fid, '--------------------------------------------------------------------');
    
    r1 = T_c(strcmp(T_c.subject, 'nagao'), :);
    r2 = T_c(strcmp(T_c.subject, 'nagao2'), :);
    
    if height(r1) == 0 || height(r2) == 0
        print_to(fid, '  One or both files missing/failed. Skipping comparison.');
        if height(r1) > 0; print_to(fid, '    nagao  format=%s', r1.format{1}); end
        if height(r2) > 0; print_to(fid, '    nagao2 format=%s', r2.format{1}); end
        continue;
    end
    
    metrics = {
        'n_dev_kept',         'HIGH', 'd';
        'rej_rate_dev',       'LOW',  'pct';
        'bl_sd_dev_uV',       'LOW',  'uV';
        'N1_Fz_uV',           'ABS',  'uV';
        'MMN_Fz_uV',          'ABS',  'uV';
        'SNR_Fz',             'HIGH', 'x';
        'split_half_topo_r',  'HIGH', 'r';
        'n_sig_ch_p05',       'HIGH', 'd';
    };
    
    print_to(fid, '');
    print_to(fid, '  format used: nagao=%s | nagao2=%s', r1.format{1}, r2.format{1});
    print_to(fid, '');
    print_to(fid, '  %-22s  %12s  %12s  %s', 'Metric', 'nagao', 'nagao2', 'winner');
    print_to(fid, '  %s', repmat('-', 1, 70));
    
    n_win_1 = 0; n_win_2 = 0; n_tie = 0;
    
    for m = 1:size(metrics, 1)
        col = metrics{m, 1};
        dir = metrics{m, 2};
        fmt = metrics{m, 3};
        v1 = r1.(col);
        v2 = r2.(col);
        
        winner = compute_winner(v1, v2, dir);
        if strcmp(winner, 'nagao');  n_win_1 = n_win_1 + 1;
        elseif strcmp(winner, 'nagao2'); n_win_2 = n_win_2 + 1;
        else; n_tie = n_tie + 1;
        end
        
        s1 = fmt_val(v1, fmt);
        s2 = fmt_val(v2, fmt);
        print_to(fid, '  %-22s  %12s  %12s  %s', col, s1, s2, winner);
    end
    
    print_to(fid, '  %s', repmat('-', 1, 70));
    print_to(fid, '  Score:                              %d wins        %d wins   (ties: %d)', ...
        n_win_1, n_win_2, n_tie);
    
    if n_win_1 > n_win_2
        cond_winner = 'nagao';
        wins.nagao = wins.nagao + 1;
    elseif n_win_2 > n_win_1
        cond_winner = 'nagao2';
        wins.nagao2 = wins.nagao2 + 1;
    else
        cond_winner = 'TIE';
        wins.tie = wins.tie + 1;
    end
    print_to(fid, '  ===> Condition %s recommended: %s', upper(cond), upper(cond_winner));
end

print_to(fid, '');
print_to(fid, '====================================================================');
print_to(fid, '   OVERALL RECOMMENDATION');
print_to(fid, '====================================================================');
print_to(fid, '  Conditions won by nagao  : %d / 4', wins.nagao);
print_to(fid, '  Conditions won by nagao2 : %d / 4', wins.nagao2);
print_to(fid, '  Ties                     : %d / 4', wins.tie);
print_to(fid, '');

if wins.nagao > wins.nagao2
    print_to(fid, '  ===> CONSIDER USING nagao  (S10 in current mapping)');
    print_to(fid, '       *** If nagao wins across most conditions, switch S10 to nagao ***');
elseif wins.nagao2 > wins.nagao
    print_to(fid, '  ===> KEEP CURRENT nagao2 as S10 (current configuration)');
else
    print_to(fid, '  ===> COIN-FLIP: tie at the metric level');
    print_to(fid, '       Inspect diagnostic PNG and decide qualitatively');
end

print_to(fid, '');
print_to(fid, 'Note: this is a metric-based recommendation. Inspect the side-by-side');
print_to(fid, 'ERP plots in nagao_quality_diagnostic.png before final decision.');
print_to(fid, '====================================================================');

fclose(fid);

end


function print_to(fid, fmt, varargin)
fprintf(fmt, varargin{:}); fprintf('\n');
fprintf(fid, fmt, varargin{:}); fprintf(fid, '\n');
end


function w = compute_winner(v1, v2, dir)
if isnan(v1) && isnan(v2); w = 'TIE-nan'; return; end
if isnan(v1); w = 'nagao2'; return; end
if isnan(v2); w = 'nagao'; return; end
switch dir
    case 'HIGH'
        if v1 > v2; w = 'nagao'; elseif v2 > v1; w = 'nagao2'; else; w = 'TIE'; end
    case 'LOW'
        if v1 < v2; w = 'nagao'; elseif v2 < v1; w = 'nagao2'; else; w = 'TIE'; end
    case 'ABS'
        if abs(v1) > abs(v2); w = 'nagao';
        elseif abs(v2) > abs(v1); w = 'nagao2';
        else; w = 'TIE';
        end
end
end


function s = fmt_val(v, fmt)
if isnan(v); s = '   NaN'; return; end
switch fmt
    case 'd',   s = sprintf('%d', round(v));
    case 'pct', s = sprintf('%.1f%%', v*100);
    case 'uV',  s = sprintf('%+.2f', v);
    case 'x',   s = sprintf('%.2f', v);
    case 'r',   s = sprintf('%.3f', v);
    otherwise,  s = sprintf('%.3g', v);
end
end


function plot_diagnostic_panel(diffwaves, ave_std_all, ave_dev_all, ...
    conditions, subjects, times_ms, ch_names, fz, cz, f4, ...
    n1_mask, mmn_mask, out_dir)

n_cond = length(conditions);
fig = figure('Position', [50 50 1800 1100], 'Name', 'nagao vs nagao2 quality');

ch_idx_list = [fz cz f4];
ch_lbls = {'Fz', 'Cz', 'F4'};
colors = {[0 0.4 0.8], [0.8 0.2 0.2]};

for ci = 1:n_cond
    for chi = 1:3
        subplot(n_cond, 3, (ci-1)*3 + chi);
        hold on;
        n1_x = times_ms(find(n1_mask, 1, 'first'));
        n1_w = times_ms(find(n1_mask, 1, 'last')) - n1_x;
        mmn_x = times_ms(find(mmn_mask, 1, 'first'));
        mmn_w = times_ms(find(mmn_mask, 1, 'last')) - mmn_x;
        
        yl = [-6 4];
        fill([n1_x n1_x+n1_w n1_x+n1_w n1_x], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.7 0.9 1], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        fill([mmn_x mmn_x+mmn_w mmn_x+mmn_w mmn_x], [yl(1) yl(1) yl(2) yl(2)], ...
            [0.7 1 0.7], 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        
        legh = [];
        leglbl = {};
        for si = 1:2
            dw = diffwaves{ci, si};
            if all(isnan(dw(:)))
                continue;
            end
            h = plot(times_ms, dw(:, ch_idx_list(chi)), ...
                'Color', colors{si}, 'LineWidth', 1.6);
            legh(end+1) = h; %#ok<AGROW>
            leglbl{end+1} = subjects{si}; %#ok<AGROW>
        end
        
        xlim([-50 350]);
        ylim(yl);
        set(gca, 'YDir', 'reverse');
        if chi == 1
            ylabel(sprintf('%s\nµV (diff)', conditions{ci}), 'FontWeight', 'bold');
        end
        if ci == n_cond
            xlabel('Time (ms)');
        end
        if ci == 1
            title(sprintf('%s   [N1=80-130 / MMN=130-200]', ch_lbls{chi}));
        else
            title(ch_lbls{chi});
        end
        if ci == 1 && chi == 3 && ~isempty(legh)
            legend(legh, leglbl, 'Location', 'northeast');
        end
        grid on; box on;
    end
end

sgtitle('nagao vs nagao2: deviant-minus-standard difference wave (4 conditions × 3 frontocentral sites)', ...
    'FontWeight', 'bold', 'FontSize', 13);

saveas(fig, fullfile(out_dir, 'nagao_quality_diagnostic.png'));
saveas(fig, fullfile(out_dir, 'nagao_quality_diagnostic.fig'));

end
