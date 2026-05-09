function compute_continuous_std_all()

input_root = uigetdir('', 'Select root folder containing all .mat files (subjects × conditions)');
if isequal(input_root, 0); return; end

[parent, ~] = fileparts(input_root);
out_input = inputdlg('Output CSV name:', 'Output', 1, {'continuous_std_n22.csv'});
if isempty(out_input); return; end
out_path = fullfile(parent, out_input{1});

reject_input = inputdlg('Bad-segment threshold (µV)?', 'Threshold', 1, {'100'});
if isempty(reject_input); return; end
reject_uv = str2double(reject_input{1});

window_input = inputdlg('Detection window (s)?', 'Window', 1, {'1.0'});
if isempty(window_input); return; end
window_s = str2double(window_input{1});

fs = 1024;
n_ch = 32;
window_n = round(window_s * fs);
[b_bp, a_bp] = butter(4, [1 20] / (fs/2), 'bandpass');

cond_suffix = {'sitting', 'free', 'water', 'nenndo'};

mat_files = dir(fullfile(input_root, '**/*.mat'));
fprintf('Found %d .mat files in %s\n\n', length(mat_files), input_root);

results = {};
for fi = 1:length(mat_files)
    fname = mat_files(fi).name;
    fpath = fullfile(mat_files(fi).folder, fname);

    cond = '';
    for ci = 1:length(cond_suffix)
        if endsWith(fname, ['_' cond_suffix{ci} '.mat'])
            cond = cond_suffix{ci};
            break;
        end
    end
    if isempty(cond)
        continue;
    end
    subj = erase(fname, ['_' cond '.mat']);

    if endsWith(subj, '_raw') || endsWith(subj, '_ch23fixed') || endsWith(subj, '_clean')
        continue;
    end

    fprintf('[%d/%d] %s (subject=%s, cond=%s)\n', fi, length(mat_files), fname, subj, cond);

    try
        d = load(fpath);
        if ~isfield(d, 'cs_new')
            fprintf('  no cs_new, skip\n');
            continue;
        end
        cs = d.cs_new;
        if size(cs, 2) < n_ch
            fprintf('  fewer than 32 ch, skip\n');
            continue;
        end

        X = cs(:, 1:n_ch) - mean(cs(:, 1:n_ch), 2);
        X = filtfilt(b_bp, a_bp, X);
        n_t = size(X, 1);

        std_raw = mean(std(X));
        max_raw = max(abs(X(:)));

        n_segs = floor(n_t / window_n);
        keep = true(n_segs, 1);
        for s = 1:n_segs
            idx = (s-1)*window_n + 1 : s*window_n;
            seg = X(idx, :);
            if max(abs(seg(:))) > reject_uv
                keep(s) = false;
            end
        end

        kept_idx = false(n_t, 1);
        for s = 1:n_segs
            if keep(s)
                kept_idx((s-1)*window_n+1 : s*window_n) = true;
            end
        end
        if sum(kept_idx) > 0
            std_clean = mean(std(X(kept_idx, :)));
        else
            std_clean = NaN;
        end
        pct_kept = sum(keep) / n_segs * 100;

        results{end+1, 1} = subj;
        results{end, 2} = cond;
        results{end, 3} = n_t / fs;
        results{end, 4} = std_raw;
        results{end, 5} = std_clean;
        results{end, 6} = max_raw;
        results{end, 7} = pct_kept;
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
    end
end

T = cell2table(results, 'VariableNames', ...
    {'subject','condition','duration_s','std_raw_uV','std_clean_uV','max_uV','pct_clean'});
writetable(T, out_path);
fprintf('\n=== Saved: %s ===\n', out_path);
fprintf('  Rows: %d\n', height(T));

fprintf('\n=== Summary by condition ===\n');
for ci = 1:length(cond_suffix)
    cond = cond_suffix{ci};
    sub = T(strcmp(T.condition, cond), :);
    fprintf('%s: n=%d, mean clean std=%.2f µV (range %.2f-%.2f)\n', ...
        cond, height(sub), mean(sub.std_clean_uV), ...
        min(sub.std_clean_uV), max(sub.std_clean_uV));
end

msgbox(sprintf(['Done. %d subject-conditions processed.\n\nSaved: %s'], ...
    height(T), out_path), 'Done');

end
