function imokawa_loo_topo_sim()

mat_dir = uigetdir('', 'imokawa .mat ファイルのフォルダ');
if isequal(mat_dir, 0); return; end
per_subject_dir = uigetdir('', 'n=21 per_subject/ フォルダ');
if isequal(per_subject_dir, 0); return; end

fs = 1024; n_pre = 51; n_post = 460;
epoch_len = n_pre + n_post + 1;
times_ms = ((0:epoch_len-1) - n_pre) / fs * 1000;
mmn_mask = times_ms >= 130 & times_ms <= 200;
bl_mask  = times_ms < 0;
reject_uv = 100;
n_ch = 32;

[b_bp, a_bp] = butter(4, [1 20]/(fs/2), 'bandpass');

cond_labels = {'sitting','free','walter','nenndo'};
cond_keys   = {'sit',    'walk','water', 'clay'  };
mat_names   = {'imokawa_sitting.mat','imokawa_free.mat', ...
               'imokawa_water.mat',  'imokawa_nenndo.mat'};

subjects_21 = {'goto','hiramoto','hirano','isoda','kenta','koga','kurayama', ...
               'murakami','nagao2','nakamura','noguchi','sasaya','suto', ...
               'tadokoro','taguchi','takashi','takeuchi','terada','uza', ...
               'yamashita','yasuyuki'};

results = struct();

for ci = 1:4
    cond   = cond_keys{ci};
    label  = cond_labels{ci};
    mfile  = fullfile(mat_dir, mat_names{ci});

    if ~exist(mfile,'file')
        fprintf('SKIP (not found): %s\n', mat_names{ci});
        continue;
    end

    % --- imokawa topography from .mat ---
    d = load(mfile);
    X = d.cs_new;
    ref = mean(X(:,1:32), 2);
    eeg = filtfilt(b_bp, a_bp, X(:,1:32) - ref);

    trig = X(:,57);
    [pk_amp, pk_idx] = findpeaks(trig, 'MinPeakDistance', 450);
    std_idx = pk_idx(pk_amp >= -161535.583 & pk_amp <= -161535.56);
    dev_idx = pk_idx(pk_amp > -161535.56);

    [ep_std, ~] = make_epochs(eeg, std_idx, n_pre, n_post, reject_uv);
    [ep_dev, ~] = make_epochs(eeg, dev_idx, n_pre, n_post, reject_uv);
    ep_std = baseline_correct(ep_std, bl_mask);
    ep_dev = baseline_correct(ep_dev, bl_mask);

    diff_imo = mean(ep_dev,3) - mean(ep_std,3);
    imo_topo = mean(diff_imo(mmn_mask,:), 1);

    % --- n=21 group mean topography ---
    group_topos = zeros(length(subjects_21), n_ch);
    n_found = 0;
    for si = 1:length(subjects_21)
        subj = subjects_21{si};
        dev_f = fullfile(per_subject_dir, sprintf('aveDev_%s_%s.csv', cond, subj));
        std_f = fullfile(per_subject_dir, sprintf('aveStd_%s_%s.csv', cond, subj));
        if ~exist(dev_f,'file') || ~exist(std_f,'file'); continue; end
        dev_ep = readmatrix(dev_f);
        std_ep = readmatrix(std_f);
        diff_s = dev_ep - std_ep;
        diff_s = diff_s - mean(diff_s(bl_mask,:),1);
        n_found = n_found + 1;
        group_topos(n_found,:) = mean(diff_s(mmn_mask,:),1);
    end
    group_topos = group_topos(1:n_found,:);
    n21_mean = mean(group_topos, 1);

    % --- LOO topo_sim ---
    r = corrcoef(imo_topo, n21_mean);
    loo_r = r(1,2);

    results.(cond).imo_topo  = imo_topo;
    results.(cond).n21_mean  = n21_mean;
    results.(cond).loo_topo_sim = loo_r;
    results.(cond).n_subjects   = n_found;

    fprintf('%s: n=21 found=%d, LOO topo_sim=%.3f\n', cond, n_found, loo_r);
end

% --- 出力 ---
conds_out = fieldnames(results);
T = table(conds_out, ...
    cellfun(@(c) results.(c).loo_topo_sim, conds_out), ...
    cellfun(@(c) results.(c).n_subjects,   conds_out), ...
    'VariableNames', {'condition','loo_topo_sim','n21_subjects'});
disp(T);

out_path = fullfile(mat_dir, 'imokawa_loo_topo_sim.csv');
writetable(T, out_path);
fprintf('保存: %s\n', out_path);
end


function [epochs, keep] = make_epochs(eeg, idx_list, n_pre, n_post, reject_uv)
n_total  = length(idx_list);
ep_len   = n_pre + n_post + 1;
n_ch     = size(eeg,2);
ep_all   = zeros(ep_len, n_ch, n_total);
keep     = false(n_total,1);
for k = 1:n_total
    s = idx_list(k) - n_pre;
    e = idx_list(k) + n_post;
    if s < 1 || e > size(eeg,1); continue; end
    ep = eeg(s:e,:);
    if any(abs(ep(:)) > reject_uv); continue; end
    ep_all(:,:,k) = ep;
    keep(k) = true;
end
epochs = ep_all(:,:,keep);
end

function out = baseline_correct(epochs, bl_mask)
out = epochs;
for k = 1:size(epochs,3)
    out(:,:,k) = epochs(:,:,k) - mean(epochs(bl_mask,:,k),1);
end
end