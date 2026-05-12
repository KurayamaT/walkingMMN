function prepare_kinematic_dataset()

motor_dir = uigetdir('', ...
    'Select motor root directory (contains free/water/nenndo subfolders)');
if isequal(motor_dir, 0); return; end

[mf, mp] = uigetfile({'*.csv', 'CSV file'}, ...
    'Select file_mapping_n22.csv', 'file_mapping_n22.csv');
if isequal(mf, 0); return; end
map_path = fullfile(mp, mf);

default_out = fullfile(fileparts(motor_dir), 'walking_MMN_dataset_public');
out_input = inputdlg('Output directory (will create kinematic/ subfolder):', ...
    'Output', 1, {default_out});
if isempty(out_input); return; end
out_dir = out_input{1};

kin_out = fullfile(out_dir, 'kinematic');
pending_out = fullfile(kin_out, '_pending_review');
if ~exist(kin_out, 'dir'); mkdir(kin_out); end
if ~exist(pending_out, 'dir'); mkdir(pending_out); end

cond_src = {'free', 'water', 'nenndo'};
cond_dst = {'walk_free', 'walk_water', 'walk_clay'};
sensors = {'wrist', 'L3'};

mapping = readtable(map_path);
n_subj = height(mapping);
n_cond = length(cond_src);
n_sens = length(sensors);
n_total = n_subj * n_cond * n_sens;

manifest = cell(n_total, 8);
row = 0;
n_ok = 0;
n_missing = 0;
n_pending = 0;

log_fid = fopen(fullfile(out_dir, 'kinematic_rename_log.txt'), 'w');
fprintf(log_fid, 'prepare_kinematic_dataset.m run log\n');
fprintf(log_fid, 'Date: %s\n', datestr(now));
fprintf(log_fid, 'Motor directory: %s\n', motor_dir);
fprintf(log_fid, 'Mapping: %s\n', map_path);
fprintf(log_fid, 'Output: %s\n\n', kin_out);

fprintf('\n========================================================\n');
fprintf(' prepare_kinematic_dataset.m\n');
fprintf('========================================================\n');
fprintf(' Motor:          %s\n', motor_dir);
fprintf(' Output:         %s\n', kin_out);
fprintf(' Pending review: %s\n', pending_out);
fprintf(' Total expected: %d files (22 subj × 3 cond × 2 sensors)\n', n_total);
fprintf('========================================================\n\n');

for si = 1:n_subj
    sid = mapping.analysis_id{si};
    orig = '';
    if ismember('original_label', mapping.Properties.VariableNames)
        orig = mapping.original_label{si};
    end
    
    for ci = 1:n_cond
        cond_folder = fullfile(motor_dir, cond_src{ci});
        if ~exist(cond_folder, 'dir')
            fprintf('Condition folder not found: %s\n', cond_folder);
            for sn = 1:n_sens
                row = row + 1;
                manifest(row, :) = {sid, orig, cond_src{ci}, cond_dst{ci}, ...
                    sensors{sn}, '', '', 'cond_folder_missing'};
                n_missing = n_missing + 1;
            end
            continue;
        end
        
        for sn = 1:n_sens
            row = row + 1;
            sensor = sensors{sn};
            
            candidates = build_candidate_paths(cond_folder, orig, sensor);
            
            found_path = '';
            for cp = 1:length(candidates)
                if exist(candidates{cp}, 'file')
                    found_path = candidates{cp};
                    break;
                end
            end
            
            dst_name = sprintf('%s_%s_%s.csv', sid, cond_dst{ci}, sensor);
            dst_path = fullfile(kin_out, dst_name);
            
            if ~isempty(found_path)
                try
                    copyfile(found_path, dst_path, 'f');
                    [~, src_base, src_ext] = fileparts(found_path);
                    src_filename = [src_base, src_ext];
                    fprintf('[%3d/%3d] %s/%s → %s  (%s)\n', ...
                        row, n_total, cond_src{ci}, src_filename, dst_name, sid);
                    fprintf(log_fid, '[OK] %s : %s → %s\n', ...
                        sid, found_path, dst_path);
                    manifest(row, :) = {sid, orig, cond_src{ci}, cond_dst{ci}, ...
                        sensor, src_filename, dst_name, 'ok'};
                    n_ok = n_ok + 1;
                catch ME
                    fprintf('[%3d/%3d] %s  COPY FAILED: %s\n', ...
                        row, n_total, sid, ME.message);
                    fprintf(log_fid, '[ERROR] %s : copy failed (%s)\n', ...
                        sid, ME.message);
                    manifest(row, :) = {sid, orig, cond_src{ci}, cond_dst{ci}, ...
                        sensor, '', '', sprintf('error:%s', ME.message)};
                    n_missing = n_missing + 1;
                end
            else
                pending_files = scan_pending_candidates(cond_folder, orig, sensor);
                
                if isempty(pending_files)
                    fprintf('[%3d/%3d] %s  %s/%s sensor=%s : MISSING (no candidates)\n', ...
                        row, n_total, sid, cond_src{ci}, orig, sensor);
                    fprintf(log_fid, '[MISSING] %s_%s_%s : no candidate file found\n', ...
                        sid, cond_dst{ci}, sensor);
                    manifest(row, :) = {sid, orig, cond_src{ci}, cond_dst{ci}, ...
                        sensor, '', '', 'missing'};
                    n_missing = n_missing + 1;
                else
                    pending_subdir = fullfile(pending_out, ...
                        sprintf('%s_%s', sid, cond_dst{ci}));
                    if ~exist(pending_subdir, 'dir'); mkdir(pending_subdir); end
                    
                    copied = {};
                    for pf = 1:length(pending_files)
                        [~, b, e] = fileparts(pending_files{pf});
                        dst_pending = fullfile(pending_subdir, [b e]);
                        try
                            copyfile(pending_files{pf}, dst_pending, 'f');
                            copied{end+1} = [b e]; %#ok<AGROW>
                        catch
                        end
                    end
                    
                    fprintf('[%3d/%3d] %s  %s/%s sensor=%s : PENDING REVIEW (%d candidates)\n', ...
                        row, n_total, sid, cond_src{ci}, orig, sensor, length(copied));
                    for cc = 1:length(copied)
                        fprintf('              candidate: %s\n', copied{cc});
                    end
                    fprintf(log_fid, '[PENDING] %s_%s_%s : %d candidates copied to %s\n', ...
                        sid, cond_dst{ci}, sensor, length(copied), pending_subdir);
                    
                    manifest(row, :) = {sid, orig, cond_src{ci}, cond_dst{ci}, ...
                        sensor, strjoin(copied, ';'), '', 'pending_review'};
                    n_pending = n_pending + 1;
                end
            end
        end
    end
end

T = cell2table(manifest, 'VariableNames', ...
    {'analysis_id', 'original_label', 'source_condition', 'output_condition', ...
    'sensor', 'source_file', 'output_file', 'status'});
manifest_path = fullfile(out_dir, 'kinematic_manifest.csv');
writetable(T, manifest_path);

dst_files = dir(fullfile(kin_out, 'S*_*_*.csv'));
n_present = length(dst_files);

fprintf('\n========================================================\n');
fprintf(' SUMMARY\n');
fprintf('========================================================\n');
fprintf(' Expected:                   %d files\n', n_total);
fprintf(' Successfully copied:        %d\n', n_ok);
fprintf(' Missing (no candidate):     %d\n', n_missing);
fprintf(' Pending review:             %d\n', n_pending);
fprintf(' Files in kinematic/:        %d\n', n_present);
fprintf(' Manifest:                   %s\n', manifest_path);
fprintf(' Pending review folder:      %s\n', pending_out);
fprintf('========================================================\n');

fprintf(log_fid, '\nSummary: %d expected, %d ok, %d missing, %d pending\n', ...
    n_total, n_ok, n_missing, n_pending);
fclose(log_fid);

if n_pending > 0
    msgbox(sprintf(['Done.\n\n' ...
        'OK:      %d\n' ...
        'Missing: %d\n' ...
        'Pending: %d  ← inspect %s\n\n' ...
        'Open each pending subfolder, judge usability, and\n' ...
        'manually copy usable files into kinematic/ with the\n' ...
        'naming pattern Sxx_walk_<cond>_<wrist|L3>.csv'], ...
        n_ok, n_missing, n_pending, pending_out), ...
        'Done (pending review needed)', 'warn');
else
    msgbox(sprintf(['All %d expected files copied successfully.\n\n' ...
        'Output: %s'], n_ok, kin_out), 'Done');
end

end


function paths = build_candidate_paths(cond_folder, orig_label, sensor)
paths = {};
if isempty(orig_label); return; end

switch sensor
    case 'wrist'
        paths{end+1} = fullfile(cond_folder, sprintf('%s.csv', orig_label));
        paths{end+1} = fullfile(cond_folder, sprintf('%s.CSV', orig_label));
    case 'L3'
        paths{end+1} = fullfile(cond_folder, sprintf('%sC.csv', orig_label));
        paths{end+1} = fullfile(cond_folder, sprintf('%s_C.csv', orig_label));
        paths{end+1} = fullfile(cond_folder, sprintf('%sc.csv', orig_label));
        paths{end+1} = fullfile(cond_folder, sprintf('%sC.CSV', orig_label));
end
end


function found = scan_pending_candidates(cond_folder, orig_label, sensor)
found = {};
if isempty(orig_label); return; end

all_csvs = dir(fullfile(cond_folder, '*.csv'));

prefixes = {orig_label};
if length(orig_label) > 0
    prefixes{end+1} = lower(orig_label);
end

for k = 1:length(all_csvs)
    fname = all_csvs(k).name;
    [~, base, ~] = fileparts(fname);
    
    matched = false;
    for p = 1:length(prefixes)
        pf = prefixes{p};
        if startsWith(lower(base), lower(pf))
            matched = true;
            break;
        end
    end
    
    if matched
        switch sensor
            case 'wrist'
                if ~endsWith(base, 'C', 'IgnoreCase', true) && ...
                   ~endsWith(base, '_C', 'IgnoreCase', true)
                    found{end+1} = fullfile(cond_folder, fname); %#ok<AGROW>
                end
            case 'L3'
                if endsWith(base, 'C', 'IgnoreCase', true) || ...
                   endsWith(base, '_C', 'IgnoreCase', true)
                    found{end+1} = fullfile(cond_folder, fname); %#ok<AGROW>
                end
        end
    end
end
end
