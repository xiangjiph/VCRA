%% Find region in Allen atlas
% Find all the substructure under a structure
% tmp_region_ind = allen_atlas.fun_name_to_ind(allen_atlas.structure_table.name, 'Piriform')
% tmp_region_id = allen_atlas.structure_table.id(tmp_region_ind);
%
% allen_atlas.structure_table(tmp_region_ind, :)
%
% allen_atlas.structure_table(any(allen_atlas.structure_id_path_mat == 987, 2) & allen_atlas.structure_id_depth == 6, :)
% Preprocess the table
% Set the regional id list according to the xlxs
% glucose_metabolism_data_fp = 'Metadata/Rodent_brain_regional_glucose_metabolism.xlsx';
% glucose_metabolism_data = readtable(glucose_metabolism_data_fp);
% glucose_metabolism_data.has_valid_atals_id_Q = true(size(glucose_metabolism_data, 1), 1);
% for iter_str = 1 : size(glucose_metabolism_data, 1)
%    tmp_id = glucose_metabolism_data.AllenAtlasID{iter_str};
%    tmp_id = strsplit(tmp_id, ',');
%    tmp_id = cellfun(@str2double, tmp_id, 'UniformOutput', false);
%    tmp_id = cat(1, tmp_id{:});
%    glucose_metabolism_data.AllenAtlasID{iter_str} = tmp_id;
%    if any(isnan(tmp_id))
%        glucose_metabolism_data.has_valid_atals_id_Q(iter_str) = false;
%    end
% end
% glucose_metabolism_data = glucose_metabolism_data(glucose_metabolism_data.has_valid_atals_id_Q, :);
% glucose_metabolism_data = table2struct(glucose_metabolism_data, 'ToScalar', true);
%
% glucose_metabolism_data.Structure = cellfun(@(x) strrep(x, 'n.', 'nucleus'), glucose_metabolism_data.Structure, 'UniformOutput', false);
% glucose_metabolism_data_fp = strrep(glucose_metabolism_data_fp, '.xlsx', '_in_atlas.mat');
% save(glucose_metabolism_data_fp, '-struct', 'glucose_metabolism_data');
%%
clc;clear;close all;
%% Parameters
brainTissueDensity = 1.042; % g/ml
% estimated from figure of DiResta et al 1991 Measurement of brain
% tissue specific gravity using pycnometry
% convert the unit of glucose metabolism to mol/(s * L)
coeffGlu2SI = 1e-6/(60 * 100 / (brainTissueDensity * 1e3)); 
glucoseOxygenIndex = 5.65;
alphaO2 = 9.8214e-4 / 760;
% Oxygen diffusion coefficient
% 2.78e-9 m^2/s in water at 40.2 C
% 2.52e-9 m^2/s in water at 35.1 C
diffCoeffO2 = 1.9e-9; % m^2/s in rodent brain (Clark et al 1978)
% Oxygen solubility
% 1.3e-3 mol(L * atm) in water;
% 9.8214e-4 mol/(L * atm) in small rodent brain (in vivo) Clark et al 1978
% alphaO2 = 1.3e-3 / 760; % mol/(L * mmHg) in water
%% Settings
de_opt_str = struct;
de_opt_str.cube_stat_Q = true;
de_opt_str.node_feature_Q = false;
de_opt_str.link_feature_Q = false;
de_opt_str.vessel_graph_Q = false;
de_opt_str.depth_dependent_density_Q = false;
de_opt_str.merge_cc_Q = true;
de_opt_str.save_Q = false;
de_opt_str.save_fp = [];

if de_opt_str.merge_cc_Q
    num_cc = 1;
else
    num_cc = 2;
end
%%  Load registration and brain mask
DataManager = FileManager;
dataset_name = 'WholeBrain';
image_grid_version = '240_cube';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
skel_version = '240_cube_rec';
merge_stack_name = 'all_stack';
pO2_folder_name = 'pO2_SA_li_sc';
reconstruction_version = '240_cube_recon_sc';
wb_stat_folder_name = 'whole_brain_stat_sc';
plot_stack_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);

num_stack = numel(stack_list);
allen_atlas = load('Allen_atlas.mat');
registration_version = 'Allen_2017_25um_landmark.mat';
Allen_atlas_id = load('Allen_atlas_id.mat');
tic_load_data = tic;

% Fitting parameters
pO2_SA_filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, ...
    merge_stack_name), 'pO2_lm_SA_result_ptl.mat');
pO2_SA = DataManager.load_data(pO2_SA_filepath);

min_in_cc_vol_f = 0.75;
min_in_brain_mask_ratio = 1;
min_cap2vsl_vol_r = pO2_SA.min_cap2vsl_vol_fraction;

upper_prctile = pO2_SA.selection_upper_prctile;
lower_prctile = pO2_SA.selection_lower_prctile;
%% Load data
wb_data_cell = cell(num_stack, 1);
for iter_stack = 1 : num_stack
    tmp_tic = tic;
    tmp_stack = stack_list{iter_stack};
    wb_region_stat_data = fun_analysis_load_whole_brain_data_for_regional_analysis(...
        dataset_name, tmp_stack, image_grid_version, reconstruction_version, ...
        skel_version, registration_version, wb_stat_folder_name, false);
    wb_region_stat_data.cube_stat = rmfield(wb_region_stat_data.cube_stat, {'wb_ai_all_lw', 'wb_ai_all_vw', 'wb_ai_cap_lw', 'wb_ai_cap_vw', 'node_stat'});
    
    tmp_pO2_stat_filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, tmp_stack), ...
        sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, tmp_stack, skel_version, pO2_folder_name));
    wb_region_stat_data.pO2_data = DataManager.load_data(tmp_pO2_stat_filepath);
    %% Add selection
    wb_region_stat_data.is_internal_cube_Q = wb_region_stat_data.cube_stat.cube_in_brain_mask_ratio >= min_in_brain_mask_ratio;
    tmp_primary_selection_Q = wb_region_stat_data.is_internal_cube_Q & ...
        wb_region_stat_data.cube_stat.cap2vsl_vol_fraction >= min_cap2vsl_vol_r;
    wb_region_stat_data.is_dominanted_by_capillary_Q = tmp_primary_selection_Q;
    
    wb_region_stat_data.is_selected_cube_Q = pO2_SA.is_selected_cube_Q{iter_stack};
%         fun_analysis_is_in_percentile_interval_Q(wb_region_stat_data.cube_stat.link_all_stat.mean.nearest_tissue_dt_max, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
%         fun_analysis_is_in_percentile_interval_Q(wb_region_stat_data.cube_stat.link_all_stat.num_data.length, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
%         fun_analysis_is_in_percentile_interval_Q(wb_region_stat_data.cube_stat.link_length_density_m_mm3, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
%         fun_analysis_is_in_percentile_interval_Q(wb_region_stat_data.pO2_data.local_dt_stat.mean, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
%         fun_analysis_is_in_percentile_interval_Q(wb_region_stat_data.pO2_data.local_pO2_stat.mean, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
%         fun_analysis_is_in_percentile_interval_Q(wb_region_stat_data.pO2_data.pO2_lm.pO2_mean(:, 11), lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
%         fun_analysis_is_in_percentile_interval_Q(wb_region_stat_data.pO2_data.dt_lm.dt_mean(:, 11), lower_prctile, upper_prctile, tmp_primary_selection_Q);
    %%
    wb_data_cell{iter_stack} = wb_region_stat_data;
    clearvars wb_region_stat_data
    fprintf('Finish loading grid, registration and brain mask distance transform for stack %s. Elapsed time is %f seconds\n', tmp_stack, toc(tmp_tic));
end
fprintf('Finish loading data. Elapsed time is %f seconds.\n', toc(tic_load_data));
% Load metabolism data
allen_atlas_map_old_id_to_new_id = wb_data_cell{1}.registration.map_oldID_to_newID;
glucose_metabolism_data = load('Metadata/Rodent_brain_regional_glucose_metabolism_in_atlas.mat');
glucose_metabolism_table = struct2table(glucose_metabolism_data);
%% Fitting parameters
local_extrema_window_size_um = wb_data_cell{iter_stack}.pO2_data.local_extrema_window_size;
num_window_size = numel(local_extrema_window_size_um);
%% Set visualization group information
vis_info = cell(3, 0);
vis_info(:, end+1) = {'Awake (Unstressed)', 'Bryan et al 1983 (Cerebral glucose utilization in awake unstressed rats)', {'Medial septum'}};
vis_info(:, end+1) = {'Awake (Unstressed)', 'Hawkins et al AM J Physiology 1985', {}};
% vis_info(:, end+1) = {'Awake', 'Sokoloff JCBFM 1981', {}};
% vis_info(:, end+1) = {'Alert', 'Hawkins et al Stroke 1979', {}};
% vis_info(:, end+1) = {'Awake (Stressed, immobilized)', 'Bryan et al 1983 (Cerebral glucose utilization in awake unstressed rats)', {}};
% vis_info(:, end+1) = {'Active', 'Jay et al Brain Research 1985', {}};
% vis_info(:, end+1) = {'Drowsy', 'Jay et al Brain Research 1985', {}};
num_reference_data = size(vis_info, 2);
% exclude_data_for_fitting = {2, 'Corpus callosum'};
exclude_data_for_fitting = [];
has_excluded_data_Q = ~isempty(exclude_data_for_fitting);
if has_excluded_data_Q
    vis_subfolder_name = sprintf('Kglu_fit_%s_eo', pO2_folder_name);
else
    vis_subfolder_name = sprintf('Kglu_fit_%s', pO2_folder_name);
end
%% Set the regional id list according to the xlxs
% Select part of the data for fitting
for iter_vis = 1 : num_reference_data
    % Select the data from table 
    select_condition = vis_info{1, iter_vis};
    select_paper = vis_info{2, iter_vis};
    exclude_region_name = vis_info{3, iter_vis};
    
    selected_data_Q = cellfun(@(x) strcmp(x, select_condition),  glucose_metabolism_data.Condition) & ...
        cellfun(@(x) strcmp(x, select_paper), glucose_metabolism_data.Reference);
    
    selected_data_ind = find(selected_data_Q);
    analysis_region_id_list = glucose_metabolism_data.AllenAtlasID(selected_data_Q);
    analysis_region_id_name_list = glucose_metabolism_data.Structure(selected_data_Q);
    selected_data_reference = unique(glucose_metabolism_data.Reference(selected_data_Q));
    filename_prefix = strsplit(selected_data_reference{:}, {' ', '_'});
    filename_prefix = sprintf('%s_%s', filename_prefix{1}, strrep(select_condition, ' ', '_'));
    
    visualization_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
        merge_stack_name), vis_subfolder_name, filename_prefix);
    num_region = numel(analysis_region_id_list);
    if ~isfolder(visualization_folder)
        mkdir(visualization_folder);
    end
    %% Collect all the region statistics in all the stack
    region_data = cell(num_cc, num_region, num_stack);
    collect_data_tic = tic;
    for iter_stack = 1 : num_stack
        for iter_region = 1 : num_region
            tmp_tic = tic;
            tmp_region_id = analysis_region_id_list{iter_region};
            tmp_region_id = full(allen_atlas_map_old_id_to_new_id(tmp_region_id));
            if ~all(tmp_region_id > 0)
                warning('Invalid atlas ID. Skip');
                continue;
            end
            tmp_region_name = analysis_region_id_name_list{iter_region};
            tmp_wb_data = wb_data_cell{iter_stack};
            tmp_region_data = fun_analysis_get_atlas_regional_data(tmp_wb_data, ...
                tmp_region_id, tmp_region_name, de_opt_str);                      
            %% Process region pO2 data
            tmp_pO2_data_field_name = setdiff(fieldnames(tmp_wb_data.pO2_data), ...
                {'dataset_name', 'stack', 'filepath', 'stack', 'skel_version', 'pO2_data_folder', 'local_extrema_window_size'});
            if all(isfield(tmp_region_data, {'cc_1', 'cc_2'})) || isfield(tmp_region_data, 'cc_merged')
                for iter_cc = 1 : tmp_region_data.num_cc
                    switch tmp_region_data.num_cc
                        case 1
                            tmp_extract_field_name = 'cc_merged';
                        case 2
                            tmp_extract_field_name = sprintf('cc_%d', iter_cc);
                    end
                    tmp_cc_data = tmp_region_data.(tmp_extract_field_name);
                    tmp_cc_data.local_cube_stat = rmfield(tmp_cc_data.local_cube_stat, {'mask_volume', 'num_surface_voxel'});
                    % Further select the data 
                    tmp_cc_data.is_selected_Q = tmp_wb_data.is_selected_cube_Q(tmp_cc_data.is_in_cc_cube_Q);                    
                    
                    for iter_field_name = 1 : numel(tmp_pO2_data_field_name)
                        tmp_field_name = tmp_pO2_data_field_name{iter_field_name};
                        tmp_cc_data.pO2_data.(tmp_field_name) = fun_structure_field_indexing(tmp_wb_data.pO2_data.(tmp_field_name), ...
                            tmp_cc_data.is_in_cc_cube_Q);
                    end
                    region_data{iter_cc, iter_region, iter_stack} = tmp_cc_data;
                end
            end
            toc(tmp_tic);
        end
    end
    fprintf('Finish collecting all the region data. Elapsed time is %f seconds.\n', toc(collect_data_tic));
    %% Compute the average capillary length density for each image stack
    region_stat = cell(num_region, num_stack);
    
    cc_merge_field_name = {'is_selected_Q', 'local_cube_stat.in_cc_vol_f', 'local_cube_stat.cube_in_brain_mask_ratio', ...
        'local_cube_stat.capillary_length_density_m_mm3', 'local_cube_stat.link_length_density_m_mm3', ...
        'local_cube_stat.cap2vsl_vol_fraction', ...
        'pO2_data.local_dt_stat.mean', 'pO2_data.local_dt_stat.median',...
        'pO2_data.local_pO2_stat.mean', 'pO2_data.local_pO2_stat.median', ...
        'pO2_data.dt_lm.dt_mean', 'pO2_data.dt_lm.dt_median', ...
        'pO2_data.dt_lm.pO2_mean', 'pO2_data.dt_lm.pO2_median', ...
        'pO2_data.pO2_lm.dt_mean', 'pO2_data.pO2_lm.dt_median',...
        'pO2_data.pO2_lm.pO2_mean', 'pO2_data.pO2_lm.pO2_median'};
    
    save_field_name = {'is_selected_Q', 'in_cc_vol_f', 'cube_in_brain_mask_ratio', ...
        'capillary_length_density_m_mm3', 'vessel_length_density_m_mm3', ...
        'cap2vsl_vol_fraction', ...
        'local_dt_mean', 'local_dt_median', ...
        'local_pO2_mean', 'local_pO2_median', ...
        'dt_lm_dt_v_mean', 'dt_lm_dt_v_median', ...
        'dt_lm_pO2_v_mean', 'dt_lm_pO2_v_median', ...
        'pO2_lm_dt_v_mean', 'pO2_lm_dt_v_median', ...
        'pO2_lm_pO2_v_mean', 'pO2_lm_pO2_v_median'};
    for iter_region = 1 : num_region
        for iter_stack = 1 : num_stack
            tmp_region_stat_str = struct;
            % Merge cc data here (if it was not merged previously)
            tmp_region_data = region_data(:, iter_region, iter_stack);
            if isempty(tmp_region_data{1})
                continue;
            end
            for iter_field = 1 : numel(cc_merge_field_name)
                tmp_field_name = cc_merge_field_name{iter_field};
                tmp_field_data = cellfun(@(x) fun_getfield(x, tmp_field_name), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_field_data = cat(1, tmp_field_data{:});                
                tmp_region_stat_str.data.(save_field_name{iter_field}) = tmp_field_data;
            end
            
            if isfield(tmp_region_stat_str.data, 'is_selected_Q') && false
                % Do not use the same selection here. Previously we
                % focused on the scaling behavior in the whole brian
                % scale and therefore want to remove outliers. The not
                % selected cubes were not necessarily outliers. E.g.
                % inferior colliculus and cochlear nucleus has high
                % vessel density, > 99.5% among all the 240-cubes
                % inside the brain.
                    
                tmp_is_valid = tmp_region_stat_str.data.is_selected_Q & ...
                    tmp_region_stat_str.data.in_cc_vol_f >= min_in_cc_vol_f;
            else
                % Previous version. Select 240-cube here
                tmp_is_valid = tmp_region_stat_str.data.in_cc_vol_f >= min_in_cc_vol_f & ...
                    tmp_region_stat_str.data.cube_in_brain_mask_ratio >= min_in_brain_mask_ratio & ...
                    tmp_region_stat_str.data.cap2vsl_vol_fraction > min_cap2vsl_vol_r;
            end
            tmp_region_stat_str.num_valid_cube = nnz(tmp_is_valid);
            if tmp_region_stat_str.num_valid_cube > 0
                for iter_field = 1 : numel(save_field_name)
                    tmp_field_name = save_field_name{iter_field};
                    tmp_field_data = tmp_region_stat_str.data.(tmp_field_name)(tmp_is_valid, :);
                    if size(tmp_field_data, 2) == 1
                        tmp_region_stat_str.stat.(tmp_field_name) = fun_analysis_get_basic_statistics(tmp_field_data, true);
                    else
                        tmp_region_stat_str.stat.(tmp_field_name) = fun_analysis_get_basic_statistics_in_column(tmp_field_data);
                    end
                end
                %% Deal with interpolation here                
                tmp_y_median = cellfun(@(x) fun_getfield(x, 'pO2_data.pO2_stat_in_dt_bin.y_median'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_y_median = cat(1, tmp_y_median{:});
                
                tmp_y_mean = cellfun(@(x) fun_getfield(x, 'pO2_data.pO2_stat_in_dt_bin.y_mean'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_y_mean = cat(1, tmp_y_mean{:});
                
                tmp_y_std = cellfun(@(x) fun_getfield(x, 'pO2_data.pO2_stat_in_dt_bin.y_std'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_y_std = cat(1, tmp_y_std{:});
                
                tmp_x_bin_val = cellfun(@(x) fun_getfield(x, 'pO2_data.pO2_stat_in_dt_bin.x_bin_val'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_x_bin_val = cat(1, tmp_x_bin_val{:});
                
                tmp_region_stat_str.pO2_vs_dt.median = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                    tmp_x_bin_val(tmp_is_valid), tmp_y_median(tmp_is_valid));
                tmp_region_stat_str.pO2_vs_dt.mean = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                    tmp_x_bin_val(tmp_is_valid), tmp_y_mean(tmp_is_valid));
                tmp_region_stat_str.pO2_vs_dt.std = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                    tmp_x_bin_val(tmp_is_valid), tmp_y_std(tmp_is_valid));
                %% pO2
                tmp_local_pO2_x = cellfun(@(x) fun_getfield(x, 'pO2_data.local_pO2_stat.hist_edge'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_local_pO2_x = cat(1, tmp_local_pO2_x{:});
                tmp_local_pO2_x = cellfun(@(x) x(2:end), tmp_local_pO2_x, 'UniformOutput', false);
                
                tmp_local_pO2_prctile = cellfun(@(x) fun_getfield(x, 'pO2_data.local_pO2_stat.hist_cdf'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_local_pO2_prctile = cat(1, tmp_local_pO2_prctile{:});
                
                tmp_region_stat_str.local_pO2_percentile = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                    tmp_local_pO2_x(tmp_is_valid), tmp_local_pO2_prctile(tmp_is_valid), [], 'nearest');
                
                % DT
                tmp_local_dt_x = cellfun(@(x) fun_getfield(x, 'pO2_data.local_dt_stat.hist_edge'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_local_dt_x = cat(1, tmp_local_dt_x{:});
                tmp_local_dt_x = cellfun(@(x) x(2:end), tmp_local_dt_x, 'UniformOutput', false);
                
                tmp_local_dt_prctile = cellfun(@(x) fun_getfield(x, 'pO2_data.local_dt_stat.hist_cdf'), ...
                    tmp_region_data, 'UniformOutput', false);
                tmp_local_dt_prctile = cat(1, tmp_local_dt_prctile{:});
                
                tmp_region_stat_str.local_dt_percentile = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                    tmp_local_dt_x(tmp_is_valid), tmp_local_dt_prctile(tmp_is_valid), [], 'nearest');
                % Inversed percentile interpolation
                tmp_region_stat_str.local_pO2_inv_cdf = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                    tmp_local_pO2_prctile(tmp_is_valid), tmp_local_pO2_x(tmp_is_valid), [], 'nearest');
                tmp_region_stat_str.local_dt_inv_cdf = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                    tmp_local_dt_prctile(tmp_is_valid), tmp_local_dt_x(tmp_is_valid), [], 'nearest');
            end
            region_stat{iter_region, iter_stack} = tmp_region_stat_str;
        end
    end
    %% Further select data
    min_num_bbox = 15;
    all_stack_valid_ind = find(all(~cellfun(@isempty, region_stat), 2));
    all_stack_valid_cell = region_stat(all_stack_valid_ind, :);
    num_data = cellfun(@(x) x.num_valid_cube, all_stack_valid_cell);
    
    has_many_data_Q = all(num_data > min_num_bbox, 2);
    % Remove Medial septum
    % extra_rm_Q = strcmp(analysis_region_id_name_list(all_stack_valid_ind), 'Inferior colliculus');
    if ~isempty(exclude_region_name)
        [~, exclude_ind] = intersect(analysis_region_id_name_list(all_stack_valid_ind), exclude_region_name);
        has_many_data_Q(exclude_ind) = false;
    end
    all_stack_valid_cell = all_stack_valid_cell(has_many_data_Q, :);
    all_stack_valid_ind = all_stack_valid_ind(has_many_data_Q);
    all_stack_valid_reigon_name = analysis_region_id_name_list(all_stack_valid_ind);
    all_stack_metabolism_table_ind = selected_data_ind(all_stack_valid_ind);
    % rm_structure_name = analysis_region_id_name_list(~has_many_data_Q);
    fprintf('Finish selecting data in valid regions\n');
    %% Save the table containing inforamtion for generating the plot
    data_table = struct;
    data_table.region_name = analysis_region_id_name_list(all_stack_valid_ind);
    if has_excluded_data_Q
        data_used_for_fitting_Q = true(numel(data_table.region_name), num_stack);
        [~, exclude_region_idx, ~] = intersect(data_table.region_name, exclude_data_for_fitting{2});
        data_used_for_fitting_Q(exclude_region_idx, exclude_data_for_fitting{1}) = false;
        data_table.included_for_fitting_Q = data_used_for_fitting_Q;
    end
    % Glucose metabolism
    data_table.k_glu_mean_M_per_sec = glucose_metabolism_data.GlucoseUtilizationMean_umol_min_100g_(all_stack_metabolism_table_ind) * coeffGlu2SI;
    data_table.k_glu_std_M_per_sec = glucose_metabolism_data.GlucoseUtilizationStd_umol_min_100g_(all_stack_metabolism_table_ind)  * coeffGlu2SI;
    
    data_table.eta_mean = data_table.k_glu_mean_M_per_sec * glucoseOxygenIndex ./ (alphaO2 * diffCoeffO2 * 1e12);
    data_table.eta_std = data_table.k_glu_std_M_per_sec * glucoseOxygenIndex ./ (alphaO2 * diffCoeffO2 * 1e12);
    data_table.eta_std_n = data_table.eta_std ./ data_table.eta_mean;
    
    % Collect regional statistics
    gather_info_cell = cell(2, 0);
    gather_info_cell(:, end+1) = {'capillary_length_density_m_mm3', 'cap_len_density_m_mm3'}; %#ok<*SAGROW>
    gather_info_cell(:, end+1) = {'vessel_length_density_m_mm3', 'vsl_len_density_m_mm3'};
    gather_info_cell(:, end+1) = {'pO2_lm_dt_v_mean', 'avg_ulm_dt'};
    gather_info_cell(:, end+1) = {'dt_lm_dt_v_mean', 'avg_dtlm'};
    gather_info_cell(:, end+1) = {'local_dt_mean', 'avg_dt'};
    gather_info_cell(:, end+1) = {'local_dt_median', 'med_dt'};
    gather_stat_field_name = {'mean', 'median', 'std', 'prctile_val'};
    num_stat_field_name = numel(gather_stat_field_name);
    for iter_cell = 1 : size(gather_info_cell, 2)
        tmp_gather_field = gather_info_cell{1, iter_cell};
        tmp_output_field_prefix = gather_info_cell{2, iter_cell};
        for iter_stat_field = 1 : num_stat_field_name
            tmp_stat_field_name = gather_stat_field_name{iter_stat_field};
            tmp_data = cellfun(@(x) fun_getfield(x, sprintf('stat.%s.%s', tmp_gather_field, tmp_stat_field_name)), ...
                all_stack_valid_cell, 'UniformOutput', false);
            
            if ~contains(tmp_stat_field_name, 'prctile')
                if isscalar(tmp_data{1})
                    tmp_data = cell2mat(tmp_data);
                    tmp_write_field_name = sprintf('%s_%s', tmp_output_field_prefix, tmp_stat_field_name);
                    data_table.(tmp_write_field_name) = tmp_data;
                elseif isrowvector(tmp_data{1})
                    tmp_data_cat = cell(1, num_stack);
                    for iter_cell_col = 1 : num_stack
                        tmp_data_cat{iter_cell_col} = cat(1, tmp_data{:, iter_cell_col});
                    end
                    for iter_col = 1 : size(tmp_data_cat{1}, 2)
                        tmp_save_data = cellfun(@(x) x(:, iter_col), tmp_data_cat, ...
                            'UniformOutput', false);
                        tmp_save_data = cat(2, tmp_save_data{:});
                        tmp_write_field_name = sprintf('%s_%d_%s', tmp_output_field_prefix, iter_col, tmp_stat_field_name);
                        data_table.(tmp_write_field_name) = tmp_save_data;
                    end
                end
            else                
                if isvector(tmp_data{1})
                    for iter_prctile = [25, 75]
                        tmp_prctile_idx = find([0 0.1000 1 2.5000 5 10 25 50 75 90 95 97.5000 99 99.9000 100] ...
                            == iter_prctile);
                        tmp_prctile_val = cellfun(@(x) x(:, tmp_prctile_idx)', tmp_data, 'UniformOutput', false);
                        tmp_prctile_val = cell2mat(tmp_prctile_val);
                        tmp_write_field_name = sprintf('%s_ptl%03d', tmp_output_field_prefix, iter_prctile);
                        data_table.(tmp_write_field_name) = tmp_prctile_val;
                    end
                elseif ismatrix(tmp_data{1})
                    % Percentile matrix # window size x # percentile
                    for iter_prctile = [25, 75]
                        tmp_prctile_idx = find([0 0.1000 1 2.5000 5 10 25 50 75 90 95 97.5000 99 99.9000 100] ...
                            == iter_prctile);
                        tmp_prctile_val = cellfun(@(x) x(:, tmp_prctile_idx)', tmp_data, 'UniformOutput', false);
                        tmp_data_cat = cell(1, num_stack);
                        for iter_cell_col = 1 : num_stack
                            tmp_data_cat{iter_cell_col} = cat(1, tmp_prctile_val{:, iter_cell_col});
                        end
                        for iter_col = 1 : size(tmp_data_cat{1}, 2)
                            tmp_save_data = cellfun(@(x) x(:, iter_col), tmp_data_cat, ...
                                'UniformOutput', false);
                            tmp_save_data = cat(2, tmp_save_data{:});
                            tmp_write_field_name = sprintf('%s_%d_ptl%03d', tmp_output_field_prefix, iter_col, iter_prctile);
                            data_table.(tmp_write_field_name) = tmp_save_data;
                        end
                    end
                end
            end
        end
    end
    % pO2 percentile
    stat_pO2_prctile = 0.10;
    u_ptl_10_mean = cellfun(@(x) x.local_pO2_inv_cdf.y_avg_interpolation(stat_pO2_prctile), ...
        all_stack_valid_cell);
    u_ptl_10_std = cellfun(@(x) x.local_pO2_inv_cdf.y_std_interpolation(stat_pO2_prctile), ...
        all_stack_valid_cell);
    u_ptl_10_std_n = u_ptl_10_std ./ u_ptl_10_mean;
    
    pO2_ptl_10_mean = bsxfun(@times, u_ptl_10_mean, data_table.eta_mean);
    pO2_ptl_10_std = sqrt(u_ptl_10_std_n .^ 2 + ...
        repmat(data_table.eta_std_n, 1, num_stack) .^ 2) .* abs(pO2_ptl_10_mean);
    u_ptl_10_median = cellfun(@(x1) median(cellfun(@(x2) x2(0.1), ...
        x1.local_pO2_inv_cdf.y_interpolation)), all_stack_valid_cell);
    u_ptl_10_ptl_025 = cellfun(@(x1) prctile(cellfun(@(x2) x2(0.1), ...
        x1.local_pO2_inv_cdf.y_interpolation), 25), all_stack_valid_cell);
    u_ptl_10_ptl_075 = cellfun(@(x1) prctile(cellfun(@(x2) x2(0.1), ...
        x1.local_pO2_inv_cdf.y_interpolation), 75), all_stack_valid_cell);
    
    data_table.regional_pO2_value_010_ptl25 = bsxfun(@times, u_ptl_10_ptl_025, data_table.eta_mean);
    data_table.regional_pO2_value_010_ptl75 = bsxfun(@times, u_ptl_10_ptl_075, data_table.eta_mean);
    data_table.regional_pO2_value_010_median = bsxfun(@times, u_ptl_10_median, data_table.eta_mean);
    
    data_table.regional_pO2_value_010_mean = pO2_ptl_10_mean;
    data_table.regional_pO2_value_010_std = pO2_ptl_10_std;
    
    data_table.num_datapoints = num_data(has_many_data_Q, :);
    data_table = struct2table(data_table);
    [~, tmp_sort_ind] = sort(data_table.k_glu_mean_M_per_sec , 'ascend');
    disp(data_table(tmp_sort_ind, :));
    table_fp = fullfile(visualization_folder, sprintf('%s_glucose_metabolism_with_pO2DT_data.csv', filename_prefix));
    writetable(data_table(tmp_sort_ind, :), table_fp);
    fprintf('Finish writing table to %s\n', table_fp);
    %% Fit Krogh model
    krogh_fit_x_label_list = cell(3, 0);
    krogh_fit_x_label_list(:, end+1) = {'avg_dtlm', 'Median average d_{lm} (\mum)', 'dlm_2_ulm'};
    krogh_fit_x_label_list(:, end+1) = {'avg_dtlm', 'Median average d_{lm} (\mum)', 'dlm_2_udlm'};
    krogh_fit_x_label_list(:, end+1) = {'avg_ulm_dt', 'Median average d_{ulm} (\mum)','dulm_2_ulm'};
    num_fit = size(krogh_fit_x_label_list, 2);
    kglu_vs_d_fit = struct;
     %% Generate figures
    for iter_fit = 1 : num_fit
        tmp_plot_stat_name = krogh_fit_x_label_list{1, iter_fit};
        tmp_plot_x_label = krogh_fit_x_label_list{2, iter_fit};
        tmp_krogh_coeff_source = krogh_fit_x_label_list{3, iter_fit};
        tmp_fit_info_cell = cell(num_window_size, 1);
        
        for iter_wz = 1 : num_window_size
            tmp_wz_um = local_extrema_window_size_um(iter_wz);
            [fig_hdl, tmp_fit_info_cell{iter_wz}] = fun_vis_kglu_vs_d(data_table, iter_wz, tmp_plot_stat_name, tmp_plot_x_label, ...
                plot_stack_list, pO2_SA.scaled_krogh_coeff.r_cap, ...
                pO2_SA.scaled_krogh_coeff.(tmp_krogh_coeff_source));
            fig_fp = fullfile(visualization_folder, tmp_krogh_coeff_source, ...
                sprintf('%s_kglu_vs_%s_fit_%s_wz_%d_um.png',...
                filename_prefix, tmp_plot_stat_name, tmp_krogh_coeff_source, tmp_wz_um));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
            delete(fig_hdl);
        end
        kglu_vs_d_fit.(tmp_krogh_coeff_source) = cat(1, tmp_fit_info_cell{:});
        %% Sensitivity analysis
        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        yyaxis(ax_hdl, 'left');
        plot(ax_hdl, local_extrema_window_size_um, [kglu_vs_d_fit.dlm_2_ulm.Estimate], ...
            'LineWidth', 1);
        ax_hdl.YLim = [0, 25];
        ax_hdl.YLabel.String = 'Average \DeltapO_2 (mmHg)';
        yyaxis(ax_hdl, 'right');
        plot(ax_hdl, local_extrema_window_size_um, [kglu_vs_d_fit.dlm_2_ulm.RSquaredAdj], ...
            'LineWidth', 1);
        ax_hdl.YLim = [0, 1];
        ax_hdl.YLabel.String = 'R^2-Adjusted';
        ax_hdl.FontSize = 14;
        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Window size (\mum)';
        ax_hdl.Title.String = strrep(tmp_krogh_coeff_source, '_2_', ' to ');
        fig_fp = fullfile(visualization_folder, tmp_krogh_coeff_source, ...
            sprintf('%s_kglu_%s_%s_fit_delta_pO2_vs_wd_sz.png',...
            filename_prefix, tmp_plot_stat_name, tmp_krogh_coeff_source));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
        delete(fig_hdl);       
    end
    %% Fit metabolism w.r.t. capillary length density through scaling
    length_density_fit_info = cell(5, 0);
    length_density_fit_info(:, end+1) = {'vsl_len_density_m_mm3', 'dulm_2_ulm', 'vsl_den_2_d', 'avg_dulm', 'Vessel \rho_l (m/mm^3)'};
    length_density_fit_info(:, end+1) = {'vsl_len_density_m_mm3', 'dlm_2_udlm', 'vsl_den_2_d', 'avg_dlm', 'Vessel \rho_l (m/mm^3)'};
    length_density_fit_info(:, end+1) = {'vsl_len_density_m_mm3', 'dlm_2_ulm', 'vsl_den_2_d', 'avg_dlm', 'Vessel \rho_l (m/mm^3)'};
    kglu_vs_rho_fit = struct;
    color_region_name = {'Parietal cortex', 'Hippocampus', 'Corpus callosum', ...
        'Inferior colliculus', 'Superior colliculus', 'Pons'};
        %% Generate figures
    for iter_len_fit = 1 : size(length_density_fit_info, 2)
        tmp_x_stat_name = length_density_fit_info{1, iter_len_fit};
        tmp_d_to_u_fn = length_density_fit_info{2, iter_len_fit};
        tmp_rho_fn = length_density_fit_info{3, iter_len_fit};
        tmp_rho_to_d_fn = length_density_fit_info{4, iter_len_fit};
        tmp_x_label = length_density_fit_info{5, iter_len_fit};
        
        tmp_fit_info_cell = cell(num_window_size, 1);
        for iter_wz = 1 : num_window_size
            tmp_wz_um = local_extrema_window_size_um(iter_wz);
            [fig_hdl, tmp_fit_info_cell{iter_wz}] = fun_vis_kglu_vs_rho(data_table, iter_wz, ...
                tmp_x_stat_name, tmp_x_label, ...
                pO2_SA.scaled_krogh_coeff.(tmp_d_to_u_fn), pO2_SA.scaled_krogh_coeff.r_cap, ...
                pO2_SA.(tmp_rho_fn).(tmp_rho_to_d_fn), color_region_name);
            fig_fp = fullfile(visualization_folder, tmp_d_to_u_fn, ...
                sprintf('%s_kglu_vs_%s_fit_through_%s_wz_%d_um.png', ...
                filename_prefix, tmp_rho_fn, tmp_d_to_u_fn, tmp_wz_um));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
            delete(fig_hdl);
        end
        %%
        tmp_fit_info_cell = cat(1, tmp_fit_info_cell{:});
        kglu_vs_rho_fit.(tmp_rho_fn).(tmp_d_to_u_fn) = tmp_fit_info_cell;
        
        fig_hdl = figure;
        ax_hdl = axes(fig_hdl);
        yyaxis(ax_hdl, 'left');
        plot(ax_hdl, local_extrema_window_size_um, [tmp_fit_info_cell.Estimate], ...
            'LineWidth', 1);
        ax_hdl.YLim = [0, 25];
        ax_hdl.YLabel.String = 'Average \DeltapO_2 (mmHg)';
        yyaxis(ax_hdl, 'right');
        plot(ax_hdl, local_extrema_window_size_um, [tmp_fit_info_cell.RSquaredAdj], ...
            'LineWidth', 1);
        ax_hdl.YLim = [0, 1];
        ax_hdl.YLabel.String = 'R^2-Adjusted';
        ax_hdl.FontSize = 14;
        grid(ax_hdl, 'on');
        ax_hdl.XLabel.String = 'Window size (\mum)';
        ax_hdl.Title.String = sprintf('%s to %s', erase(tmp_x_label, '(m/mm^3)'), ...
            strrep(tmp_d_to_u_fn, '_2_', ' to '));
        fig_fp = fullfile(visualization_folder, tmp_d_to_u_fn, ...
            sprintf('%s_kglu_%s_%s_fit_delta_pO2_vs_wd_sz.png',...
            filename_prefix, tmp_rho_fn, tmp_d_to_u_fn));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
        delete(fig_hdl);      
    end    
    krogh_fit_result = struct;
    krogh_fit_result.kglu_vs_d = kglu_vs_d_fit;
    krogh_fit_result.kglu_vs_rho = kglu_vs_rho_fit;
    krogh_fit_result.filepath = fullfile(visualization_folder, ...
        sprintf('%s_krogh_fit_result.mat', filename_prefix));
    DataManager.write_data(krogh_fit_result.filepath, krogh_fit_result);
end