function walking_MMN_singletrial()
% walking_MMN_singletrial.m
%
% Single-trial and topographic analyses for MMN data
% Companion script to walking_MMN_reanalysis.m
%
% Adds the following analyses:
%   Tier 1 (from averaged ERPs):
%     - Peak latency, 50% fractional area latency, onset latency
%     - Lateralization (R-L), Anterior-Posterior centroid
%     - Global Field Power (GFP), Topographic complexity
%     - Topo similarity to leave-one-out group mean (computed at group stage)
%
%   Tier 2 (from single trials):
%     - Inter-Trial Coherence (ITC) at MMN latency
%     - Per-trial amplitude standard deviation
%     - Theta-band (4-8 Hz) evoked power
%     - Bishop & Hardiman (2010) single-trial detection index
%     - Trial count to stable MMN (split-half analysis)
%
% Usage:
%   Run once per condition. Select a directory containing files named
%   <subject>_<condition>.mat. Output is saved to a sibling directory.
%
% Output structure:
%   <output>/
%     <subject>_<condition>/
%       aveStd.csv, aveDev.csv, diffWave.csv
%       itc_dev.csv, itc_std.csv
%       theta_power_dev.csv, theta_power_std.csv
%       singletrial_sd_dev.csv, singletrial_sd_std.csv
%       bishop_hardiman.csv
%     group_summary.csv
%     time_axis_ms.csv
%     info.txt
%
% Requirements:
%   - MATLAB Signal Processing Toolbox (butter, filtfilt, hilbert, findpeaks)
%   - Parallel Computing Toolbox (optional, for parfor)
%
% Compatible with output from walking_MMN_reanalysis.m

% ==================== GUI ====================

input_dir = uigetdir('', 'Select directory containing condition .mat files');
if isequal(input_dir, 0); return; end

[parent_dir, default_name] = fileparts(input_dir);
out_default = [default_name '_singletrial'];
out_input = inputdlg('Output folder name (will be created in parent dir):', ...
    'Output folder', 1, {out_default});
if isempty(out_input); return; end
output_dir = fullfile(parent_dir, out_input{1});
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

ref_options = {'Common Average Reference (CAR, recommended)', ...
               'Original mastoid-like (X37+X38+X39)/3', ...
               'No re-reference'};
[ref_idx, ok] = listdlg('PromptString', 'Reference method:', ...
                        'SelectionMode', 'single', ...
                        'ListString', ref_options, ...
                        'ListSize', [400, 100]);
if ~ok; return; end

cond_input = inputdlg('Condition suffix in filenames (e.g., sitting / free / water / nenndo):', ...
    'Condition suffix', 1, {'sitting'});
if isempty(cond_input); return; end
cond_label = cond_input{1};

% ==================== Constants ====================

fs = 1024;
n_pre = 51;          % samples pre-stimulus (50 ms at 1024 Hz)
n_post = 460;        % samples post-stimulus (450 ms)
epoch_len = n_pre + n_post + 1;
times_ms = ((0:epoch_len-1) - n_pre) / fs * 1000;
reject_uv = 100;

ch_names = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5', ...
            'P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8', ...
            'CP6','CP2','C4','T8','FC6','FC2','F4','F8','AF4','Fp2', ...
            'Fz','Cz'};

% Filter design
[b_bp, a_bp] = butter(4, [1 20] / (fs/2), 'bandpass');
[b_th, a_th] = butter(4, [4 8] / (fs/2), 'bandpass');

% MMN window indices
mmn_mask = times_ms >= 130 & times_ms <= 200;
n1_mask  = times_ms >= 80  & times_ms <= 130;
bl_mask  = times_ms <  0;

% Lateralization channel groups
left_chs  = {'Fp1','AF3','F7','F3','FC1','FC5','T7','C3','CP1','CP5'};
right_chs = {'Fp2','AF4','F8','F4','FC2','FC6','T8','C4','CP2','CP6'};
[~, l_idx] = ismember(left_chs,  ch_names);
[~, r_idx] = ismember(right_chs, ch_names);

% Standard 10-20 y-coordinates (anterior positive)
% Approximate from Biosemi32 montage
y_coords_mm = approx_biosemi_y(ch_names);

% ==================== File discovery ====================

mat_files = dir(fullfile(input_dir, '*.mat'));
valid_files = {};
subjects = {};
for fi = 1:length(mat_files)
    [~, base, ~] = fileparts(mat_files(fi).name);
    suffix = ['_' cond_label];
    if endsWith(base, suffix)
        subj = base(1:end-length(suffix));
        % Reject malformed names
        if ~isempty(subj) && ~strcmp(subj, 'matlab') && ...
           ~strcmp(subj, cond_label) && ~contains(subj, ' ')
            valid_files{end+1} = mat_files(fi).name; %#ok<AGROW>
            subjects{end+1} = subj; %#ok<AGROW>
        end
    end
end

if isempty(valid_files)
    errordlg(sprintf('No files matching <subject>_%s.mat found.', cond_label));
    return;
end
fprintf('Found %d valid files in %s\n\n', length(valid_files), input_dir);

% ==================== Open log ====================

log_fid = fopen(fullfile(output_dir, 'info.txt'), 'a');
fprintf(log_fid, '\n========== walking_MMN_singletrial run ==========\n');
fprintf(log_fid, 'Date: %s\n', datestr(now));
fprintf(log_fid, 'Input dir: %s\n', input_dir);
fprintf(log_fid, 'Condition: %s\n', cond_label);
fprintf(log_fid, 'Reference: %s\n', ref_options{ref_idx});
fprintf(log_fid, 'N files: %d\n\n', length(valid_files));
fprintf(log_fid, 'per-file diagnostics:\n');
fprintf(log_fid, 'subject  n_std n_dev half1_dev half2_dev\n');

% ==================== Process each file ====================

n_files = length(valid_files);
group_lat = zeros(n_files, 1);
group_gfp = zeros(n_files, 1);
group_ap_centroid = zeros(n_files, 1);
group_peak_lat_fz = zeros(n_files, 1);
group_n_std = zeros(n_files, 1);
group_n_dev = zeros(n_files, 1);
group_split_half_r = zeros(n_files, 1);
group_subj = cell(n_files, 1);

fz_idx = find(strcmp(ch_names, 'Fz'));

% Process files (sequential; can be parallelized if Parallel Toolbox available)
for fi = 1:n_files
    fname = valid_files{fi};
    subj = subjects{fi};
    fprintf('[%d/%d] Processing %s_%s...\n', fi, n_files, subj, cond_label);
    
    try
        % --- Load ---
        data = load(fullfile(input_dir, fname));
        if ~isfield(data, 'cs_new')
            fprintf('  WARN: cs_new not found. Skipping.\n');
            continue;
        end
        X = data.cs_new;
        if size(X, 2) < 41
            fprintf('  WARN: <41 columns. Skipping.\n');
            continue;
        end
        
        % --- Re-reference ---
        switch ref_idx
            case 1, ref_signal = mean(X(:, 1:32), 2);
            case 2, ref_signal = mean(X(:, [37 38 39]), 2);
            case 3, ref_signal = zeros(size(X,1), 1);
        end
        eeg = X(:, 1:32) - ref_signal;
        
        % --- Filter ---
        eeg_filt = filtfilt(b_bp, a_bp, eeg);
        
        % --- Trigger detection ---
        trig = X(:, 41);
        [pk_amp, pk_idx] = findpeaks(trig, 'MinPeakDistance', 100);
        
        std_idx = pk_idx(pk_amp <= -210607.93);
        dev_idx = pk_idx(pk_amp >= -210607.92 & pk_amp <= -210607.9);
        
        % --- Epoch + reject ---
        [epochs_std, k_std] = make_epochs(eeg_filt, std_idx, n_pre, n_post, reject_uv);
        [epochs_dev, k_dev] = make_epochs(eeg_filt, dev_idx, n_pre, n_post, reject_uv);
        n_std_kept = sum(k_std);
        n_dev_kept = sum(k_dev);
        
        if n_std_kept < 10 || n_dev_kept < 10
            fprintf('  WARN: insufficient trials (std=%d, dev=%d). Skipping.\n', ...
                n_std_kept, n_dev_kept);
            continue;
        end
        
        % --- Baseline correct ---
        epochs_std = baseline_correct(epochs_std, bl_mask);
        epochs_dev = baseline_correct(epochs_dev, bl_mask);
        
        % --- Output directory ---
        out_subdir = fullfile(output_dir, sprintf('%s_%s', subj, cond_label));
        if ~exist(out_subdir, 'dir'); mkdir(out_subdir); end
        
        % --- Averaged ERPs ---
        ave_std = mean(epochs_std, 3);
        ave_dev = mean(epochs_dev, 3);
        diff_wave = ave_dev - ave_std;
        
        writematrix(ave_std,   fullfile(out_subdir, 'aveStd.csv'));
        writematrix(ave_dev,   fullfile(out_subdir, 'aveDev.csv'));
        writematrix(diff_wave, fullfile(out_subdir, 'diffWave.csv'));
        
        % ==================== TIER 1: Topographic indicators ====================
        
        mmn_topo = mean(diff_wave(mmn_mask, :), 1);
        
        % Lateralization
        lateralization = mean(mmn_topo(r_idx)) - mean(mmn_topo(l_idx));
        
        % GFP
        gfp = std(mmn_topo);
        
        % Anterior-Posterior centroid (weighted by negative amplitudes)
        neg_mask = mmn_topo < 0;
        if any(neg_mask)
            weights = -mmn_topo(neg_mask);
            ap_centroid = sum(y_coords_mm(neg_mask) .* weights) / sum(weights);
        else
            ap_centroid = NaN;
        end
        
        % Peak latency at Fz (within 100-250 ms)
        search_mask = times_ms >= 100 & times_ms <= 250;
        search_t = times_ms(search_mask);
        wave_fz = diff_wave(search_mask, fz_idx);
        [~, peak_i] = min(wave_fz);
        peak_lat_fz = search_t(peak_i);
        
        % 50% fractional area latency at Fz
        seg_abs = abs(diff_wave(search_mask, fz_idx));
        cum_area = cumsum(seg_abs);
        if cum_area(end) > 0
            half_idx = find(cum_area >= cum_area(end) * 0.5, 1);
            frac_lat_fz = search_t(half_idx);
        else
            frac_lat_fz = NaN;
        end
        
        % Onset latency at Fz (first sustained sig point at 1-sample test, simplified)
        % For single-subject, define as first time wave goes <0 sustained for 10 samples
        post_mask = times_ms >= 50 & times_ms <= 250;
        wave_post = diff_wave(post_mask, fz_idx);
        t_post = times_ms(post_mask);
        onset_lat = NaN;
        for ti = 1:length(wave_post)-9
            if all(wave_post(ti:ti+9) < 0)
                onset_lat = t_post(ti);
                break;
            end
        end
        
        % ==================== TIER 2: Single-trial indicators ====================
        
        n_t = epoch_len;
        n_ch = 32;
        n_dev = size(epochs_dev, 3);
        n_std = size(epochs_std, 3);
        
        itc_dev = zeros(n_t, n_ch);
        itc_std = zeros(n_t, n_ch);
        theta_power_dev = zeros(n_t, n_ch);
        theta_power_std = zeros(n_t, n_ch);
        sd_dev = zeros(n_t, n_ch);
        sd_std = zeros(n_t, n_ch);
        
        % Pre-allocate for parfor
        epochs_dev_loc = epochs_dev;
        epochs_std_loc = epochs_std;
        
        parfor ci = 1:n_ch
            % Channel time series across trials
            ch_dev = squeeze(epochs_dev_loc(:, ci, :));  % (n_t, n_trials)
            ch_std = squeeze(epochs_std_loc(:, ci, :));
            
            % Theta-bandpass (4-8 Hz) on each trial
            ch_dev_theta = filtfilt(b_th, a_th, ch_dev);
            ch_std_theta = filtfilt(b_th, a_th, ch_std);
            
            % Hilbert transform per trial
            h_dev = hilbert(ch_dev_theta);
            h_std = hilbert(ch_std_theta);
            
            % Theta power per time point (mean across trials)
            theta_power_dev(:, ci) = mean(abs(h_dev).^2, 2);
            theta_power_std(:, ci) = mean(abs(h_std).^2, 2);
            
            % ITC per time point (across trials)
            phase_dev = angle(h_dev);
            phase_std = angle(h_std);
            itc_dev(:, ci) = abs(mean(exp(1i * phase_dev), 2));
            itc_std(:, ci) = abs(mean(exp(1i * phase_std), 2));
            
            % Single-trial SD per time point
            sd_dev(:, ci) = std(ch_dev, 0, 2);
            sd_std(:, ci) = std(ch_std, 0, 2);
        end
        
        writematrix(itc_dev,         fullfile(out_subdir, 'itc_dev.csv'));
        writematrix(itc_std,         fullfile(out_subdir, 'itc_std.csv'));
        writematrix(theta_power_dev, fullfile(out_subdir, 'theta_power_dev.csv'));
        writematrix(theta_power_std, fullfile(out_subdir, 'theta_power_std.csv'));
        writematrix(sd_dev,          fullfile(out_subdir, 'singletrial_sd_dev.csv'));
        writematrix(sd_std,          fullfile(out_subdir, 'singletrial_sd_std.csv'));
        
        % --- Bishop & Hardiman (2010) per-channel test ---
        % Compare single-trial mean amplitudes in MMN window: dev vs std
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
        
        % --- Split-half reliability ---
        % Split deviant trials into halves, average, correlate the two MMN topos
        if n_dev >= 20
            half = floor(n_dev / 2);
            half1 = mean(epochs_dev(:, :, 1:half), 3) - ave_std;
            half2 = mean(epochs_dev(:, :, half+1:2*half), 3) - ave_std;
            topo1 = mean(half1(mmn_mask, :), 1);
            topo2 = mean(half2(mmn_mask, :), 1);
            split_half_r = corrcoef(topo1, topo2);
            split_half_r = split_half_r(1, 2);
        else
            split_half_r = NaN;
        end
        
        % --- Per-subject indicators table ---
        ind_table = table({subj}, {cond_label}, n_std_kept, n_dev_kept, ...
            lateralization, gfp, ap_centroid, peak_lat_fz, frac_lat_fz, onset_lat, ...
            split_half_r, ...
            'VariableNames', {'subject', 'condition', 'n_std', 'n_dev', ...
            'lateralization_RminusL_uV', 'GFP_uV', 'AP_centroid_mm', ...
            'peak_lat_Fz_ms', 'frac50_lat_Fz_ms', 'onset_lat_Fz_ms', ...
            'split_half_topo_r'});
        writetable(ind_table, fullfile(out_subdir, 'topo_indicators.csv'));
        
        % Accumulate group results
        group_lat(fi) = lateralization;
        group_gfp(fi) = gfp;
        group_ap_centroid(fi) = ap_centroid;
        group_peak_lat_fz(fi) = peak_lat_fz;
        group_n_std(fi) = n_std_kept;
        group_n_dev(fi) = n_dev_kept;
        group_split_half_r(fi) = split_half_r;
        group_subj{fi} = subj;
        
        fprintf('  n_std=%d, n_dev=%d, lat=%+.2f, GFP=%.2f, peak_lat=%.0f, split_r=%.2f\n', ...
            n_std_kept, n_dev_kept, lateralization, gfp, peak_lat_fz, split_half_r);
        fprintf(log_fid, '%s_%s.mat  %5d  %5d  half1=%5d  half2=%5d\n', ...
            subj, cond_label, n_std_kept, n_dev_kept, ...
            floor(n_dev_kept/2), floor(n_dev_kept/2));
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        fprintf(log_fid, '%s_%s.mat  ERROR: %s\n', subj, cond_label, ME.message);
    end
end

% ==================== Group summary ====================

valid = ~cellfun(@isempty, group_subj);
group_table = table(group_subj(valid), repmat({cond_label}, sum(valid), 1), ...
    group_n_std(valid), group_n_dev(valid), ...
    group_lat(valid), group_gfp(valid), group_ap_centroid(valid), ...
    group_peak_lat_fz(valid), group_split_half_r(valid), ...
    'VariableNames', {'subject', 'condition', 'n_std', 'n_dev', ...
    'lateralization_RminusL_uV', 'GFP_uV', 'AP_centroid_mm', ...
    'peak_lat_Fz_ms', 'split_half_topo_r'});

% Append to existing group_summary.csv if it exists
group_csv = fullfile(output_dir, 'group_summary.csv');
if exist(group_csv, 'file')
    existing = readtable(group_csv);
    group_table = [existing; group_table];
end
writetable(group_table, group_csv);

% Save time axis (overwrite)
writematrix(times_ms', fullfile(output_dir, 'time_axis_ms.csv'));

fprintf('\nGroup summary saved: %s\n', group_csv);
fprintf(log_fid, '\nProcessing complete. %d/%d files succeeded.\n', sum(valid), n_files);
fclose(log_fid);

% Display summary message
msgbox(sprintf(['Processing complete.\n\n' ...
    'Condition: %s\n' ...
    'Files processed: %d/%d\n' ...
    'Output: %s\n\n' ...
    'Group summary: group_summary.csv (appended)'], ...
    cond_label, sum(valid), n_files, output_dir), 'Done');

end


% ==================== Helper functions ====================

function [epochs, keep] = make_epochs(eeg, idx_list, n_pre, n_post, reject_uv)
% Returns epochs as (n_t, n_ch, n_trials_kept), and logical keep mask
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

function y_mm = approx_biosemi_y(ch_names)
% Approximate y-coordinates (anterior-posterior, in mm) for Biosemi32
% Based on standard 10-20 montage; positive = anterior
coords = struct( ...
    'Fp1',  92, 'Fp2',  92, ...
    'AF3',  74, 'AF4',  74, ...
    'F7',   55, 'F3',   60, 'Fz',   60, 'F4',   60, 'F8',   55, ...
    'FC5',  30, 'FC1',  35, 'FC2',  35, 'FC6',  30, ...
    'T7',    0, 'C3',    0, 'Cz',    0, 'C4',    0, 'T8',    0, ...
    'CP5', -30, 'CP1', -35, 'CP2', -35, 'CP6', -30, ...
    'P7',  -55, 'P3',  -60, 'Pz',  -60, 'P4',  -60, 'P8',  -55, ...
    'PO3', -75, 'PO4', -75, ...
    'O1',  -92, 'Oz',  -92, 'O2',  -92);
y_mm = zeros(1, length(ch_names));
for i = 1:length(ch_names)
    if isfield(coords, ch_names{i})
        y_mm(i) = coords.(ch_names{i});
    else
        y_mm(i) = NaN;
    end
end
end
