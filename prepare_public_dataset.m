function prepare_public_dataset()

base_dir = uigetdir('', ...
    'Select archive ROOT directory (contains sitting_mat/free_mat/water_mat/nenndo_mat/)');
if isequal(base_dir, 0); return; end

map_path = fullfile(base_dir, 'file_mapping_n22.csv');
if ~exist(map_path, 'file')
    [mf, mp] = uigetfile({'*.csv', 'CSV file'; '*.*', 'All files'}, ...
        'Select file_mapping_n22.csv');
    if isequal(mf, 0); return; end
    map_path = fullfile(mp, mf);
end

default_out = fullfile(fileparts(base_dir), 'walking_MMN_dataset_public');
out_input = inputdlg('Output directory (will create EEG/ subfolder):', ...
    'Output', 1, {default_out});
if isempty(out_input); return; end
out_dir = out_input{1};

if exist(out_dir, 'dir')
    answer = questdlg(sprintf(['Output directory already exists:\n%s\n\n' ...
        'Continue and overwrite existing files?'], out_dir), ...
        'Overwrite?', 'Yes', 'No', 'No');
    if ~strcmp(answer, 'Yes'); return; end
else
    mkdir(out_dir);
end

eeg_out = fullfile(out_dir, 'EEG');
if ~exist(eeg_out, 'dir'); mkdir(eeg_out); end

cond_src     = {'sitting', 'free', 'water', 'nenndo'};
cond_dst     = {'sit', 'walk_free', 'walk_water', 'walk_clay'};
cond_folders = {'sitting_mat', 'free_mat', 'water_mat', 'nenndo_mat'};
cond_cols    = {'sit_file', 'free_file', 'water_file', 'clay_file'};

mapping = readtable(map_path);
if ~ismember('analysis_id', mapping.Properties.VariableNames)
    errordlg('file_mapping_n22.csv must contain column "analysis_id"');
    return;
end
for k = 1:length(cond_cols)
    if ~ismember(cond_cols{k}, mapping.Properties.VariableNames)
        errordlg(sprintf('file_mapping_n22.csv must contain column "%s"', cond_cols{k}));
        return;
    end
end

n_subj = height(mapping);
n_cond = length(cond_src);
n_total = n_subj * n_cond;

manifest = cell(n_total, 7);
row = 0;
n_ok = 0;
n_missing = 0;

log_path = fullfile(out_dir, 'rename_log.txt');
log_fid = fopen(log_path, 'w');
fprintf(log_fid, 'prepare_public_dataset.m run log\n');
fprintf(log_fid, 'Date: %s\n', datestr(now));
fprintf(log_fid, 'Base directory: %s\n', base_dir);
fprintf(log_fid, 'Mapping CSV: %s\n', map_path);
fprintf(log_fid, 'Output: %s\n\n', eeg_out);

fprintf('\n========================================================\n');
fprintf(' prepare_public_dataset.m\n');
fprintf('========================================================\n');
fprintf(' Base:    %s\n', base_dir);
fprintf(' Mapping: %s\n', map_path);
fprintf(' Output:  %s\n', eeg_out);
fprintf(' Total to copy: %d files (%d subjects × %d conditions)\n', ...
    n_total, n_subj, n_cond);
fprintf('========================================================\n\n');

for si = 1:n_subj
    sid = mapping.analysis_id{si};
    orig_label = '';
    if ismember('original_label', mapping.Properties.VariableNames)
        orig_label = mapping.original_label{si};
    end
    
    for ci = 1:n_cond
        row = row + 1;
        src_field = mapping.(cond_cols{ci}){si};
        
        if isempty(src_field) || all(isspace(src_field))
            fprintf('[%3d/%3d] %s  cond=%s  MAPPING EMPTY\n', ...
                row, n_total, sid, cond_dst{ci});
            fprintf(log_fid, '[MISSING] %s_%s : mapping field empty\n', ...
                sid, cond_dst{ci});
            manifest{row, 1} = sid;
            manifest{row, 2} = orig_label;
            manifest{row, 3} = cond_src{ci};
            manifest{row, 4} = cond_dst{ci};
            manifest{row, 5} = src_field;
            manifest{row, 6} = '';
            manifest{row, 7} = 'mapping_empty';
            n_missing = n_missing + 1;
            continue;
        end
        
        if contains(src_field, '/') || contains(src_field, '\')
            src_path = fullfile(base_dir, src_field);
        else
            src_path = fullfile(base_dir, cond_folders{ci}, src_field);
        end
        
        if ~exist(src_path, 'file')
            alt_paths = {
                fullfile(base_dir, cond_folders{ci}, src_field), ...
                fullfile(base_dir, src_field), ...
                fullfile(base_dir, cond_folders{ci}, sprintf('%s_%s.mat', orig_label, cond_src{ci}))
            };
            found = false;
            for ap = 1:length(alt_paths)
                if exist(alt_paths{ap}, 'file')
                    src_path = alt_paths{ap};
                    found = true;
                    break;
                end
            end
            if ~found
                fprintf('[%3d/%3d] %s  cond=%s  SOURCE NOT FOUND: %s\n', ...
                    row, n_total, sid, cond_dst{ci}, src_path);
                fprintf(log_fid, '[MISSING] %s_%s : source file not found (%s)\n', ...
                    sid, cond_dst{ci}, src_path);
                manifest{row, 1} = sid;
                manifest{row, 2} = orig_label;
                manifest{row, 3} = cond_src{ci};
                manifest{row, 4} = cond_dst{ci};
                manifest{row, 5} = src_field;
                manifest{row, 6} = '';
                manifest{row, 7} = 'source_not_found';
                n_missing = n_missing + 1;
                continue;
            end
        end
        
        dst_name = sprintf('%s_%s.mat', sid, cond_dst{ci});
        dst_path = fullfile(eeg_out, dst_name);
        
        try
            copyfile(src_path, dst_path, 'f');
            fprintf('[%3d/%3d] %s → %s  (%s)\n', ...
                row, n_total, src_field, dst_name, sid);
            fprintf(log_fid, '[OK] %s_%s : %s → %s\n', ...
                sid, cond_dst{ci}, src_path, dst_path);
            manifest{row, 1} = sid;
            manifest{row, 2} = orig_label;
            manifest{row, 3} = cond_src{ci};
            manifest{row, 4} = cond_dst{ci};
            manifest{row, 5} = src_field;
            manifest{row, 6} = dst_name;
            manifest{row, 7} = 'ok';
            n_ok = n_ok + 1;
        catch ME
            fprintf('[%3d/%3d] %s  cond=%s  COPY FAILED: %s\n', ...
                row, n_total, sid, cond_dst{ci}, ME.message);
            fprintf(log_fid, '[ERROR] %s_%s : copy failed (%s)\n', ...
                sid, cond_dst{ci}, ME.message);
            manifest{row, 1} = sid;
            manifest{row, 2} = orig_label;
            manifest{row, 3} = cond_src{ci};
            manifest{row, 4} = cond_dst{ci};
            manifest{row, 5} = src_field;
            manifest{row, 6} = '';
            manifest{row, 7} = sprintf('error:%s', ME.message);
            n_missing = n_missing + 1;
        end
    end
end

T = cell2table(manifest, 'VariableNames', ...
    {'analysis_id', 'original_label', 'source_condition', 'output_condition', ...
    'source_file', 'output_file', 'status'});
manifest_path = fullfile(out_dir, 'rename_manifest.csv');
writetable(T, manifest_path);

dst_files = dir(fullfile(eeg_out, 'S*_*.mat'));
n_present = length(dst_files);

fprintf('\n========================================================\n');
fprintf(' SUMMARY\n');
fprintf('========================================================\n');
fprintf(' Expected:                   %d files\n', n_total);
fprintf(' Successfully copied:        %d\n', n_ok);
fprintf(' Missing / failed:           %d\n', n_missing);
fprintf(' Files in output directory:  %d\n', n_present);
fprintf(' Manifest CSV:               %s\n', manifest_path);
fprintf(' Log file:                   %s\n', log_path);
fprintf('========================================================\n');

fprintf(log_fid, '\n========================================\n');
fprintf(log_fid, 'Summary: %d expected, %d copied, %d missing\n', ...
    n_total, n_ok, n_missing);
fclose(log_fid);

if n_missing > 0
    msgbox(sprintf(['Done with WARNINGS.\n\n' ...
        'Expected: %d\n' ...
        'Copied:   %d\n' ...
        'Missing:  %d\n\n' ...
        'Check rename_manifest.csv for details.\n\n' ...
        'Output: %s'], ...
        n_total, n_ok, n_missing, eeg_out), 'Done with warnings', 'warn');
else
    msgbox(sprintf(['All %d files copied successfully.\n\n' ...
        'Output: %s\n\n' ...
        'Next step: re-run analysis pipeline on the renamed dataset.'], ...
        n_total, eeg_out), 'Done');
end

end
