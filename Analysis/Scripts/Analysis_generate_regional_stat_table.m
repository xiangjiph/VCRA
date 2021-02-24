clc;clear;close all;
%%
DataManager = FileManager;
dataset_name = 'WholeBrain';
image_grid_version = '240_cube';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
skel_version = '240_cube_rec';
pO2_folder_name = 'pO2_SA_li_sc';
reconstruction_version = '240_cube_recon_sc';
wb_stat_folder_name = 'whole_brain_stat_sc';

plot_stack_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
linear_scaling_factor = 1.0521;
num_stack = numel(stack_list);
allen_atlas = load('Allen_atlas.mat');
registration_version = 'Allen_2017_25um_landmark.mat';
Allen_atlas_id = load('Allen_atlas_id.mat');
tic_load_data = tic;
wb_data_cell = cell(num_stack, 1);
for iter_stack = 1 : num_stack
    tmp_tic = tic;
    tmp_stack = stack_list{iter_stack};
    wb_region_stat_data = fun_analysis_load_whole_brain_data_for_regional_analysis(...
        dataset_name, tmp_stack, image_grid_version, reconstruction_version, ...
        skel_version, registration_version, wb_stat_folder_name, false);
    wb_region_stat_data.cube_stat.num_link = wb_region_stat_data.cube_stat.link_all_stat.num_data.length;
    wb_region_stat_data.cube_stat.num_cap = wb_region_stat_data.cube_stat.link_cap_stat.num_data.length;
    
    tmp_pO2_stat_filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, tmp_stack), ...
        sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, tmp_stack, skel_version, pO2_folder_name));
    wb_region_stat_data.cube_stat.pO2_data = DataManager.load_data(tmp_pO2_stat_filepath);
    %%
    wb_data_cell{iter_stack} = wb_region_stat_data;
    clearvars wb_region_stat_data
    fprintf('Finish loading grid, registration and brain mask distance transform for stack %s. Elapsed time is %f seconds\n', tmp_stack, toc(tmp_tic));
end
fprintf('Finish loading data. Elapsed time is %f seconds.\n', toc(tic_load_data));
%% Map acrony to atlas ID
registration_str = wb_data_cell{1}.registration;
tmp_is_valid_Q = isfinite(registration_str.structure_table.id);
tmp_abbrv_to_id =  containers.Map(registration_str.structure_table.acronym(tmp_is_valid_Q), ...
    registration_str.structure_table.id(tmp_is_valid_Q));

allen_atlas_map_old_id_to_new_id = wb_data_cell{1}.registration.map_oldID_to_newID;

table_data = registration_str.structure_table;
num_region = numel(table_data.acronym);

% Add field to the data table
% Need to replace NaN with -1
str_path_id = registration_str.structure_id_path_mat;
str_path_id(isnan(str_path_id)) = -1;
[str_path_id_sort, tmp_sort_ind] = sortrows(str_path_id);
str_path_abbrv_sort = cell(size(str_path_id_sort));
for iter_elem = 1 : numel(str_path_id_sort)
    tmp_id = str_path_id_sort(iter_elem);
    if tmp_id > 0
        str_path_abbrv_sort{iter_elem} = registration_str.label_2_acronym(tmp_id);
    end
end
str_path_abbrv_sort_table = cell2table(str_path_abbrv_sort);
table_data = table_data(tmp_sort_ind, :);
table_data = cat(2, str_path_abbrv_sort, table_data);

% region_data.id = deal(nan(num_region, 1));
% region_data.Structure_name = cell(num_region, 1);
% for iter_region = 1 : num_region
%     tmp_region_abbrv = region_data.Acronym{iter_region};
%     tmp_region_abbrv = erase(tmp_region_abbrv, ' ');
%     if tmp_abbrv_to_id.isKey(tmp_region_abbrv)
%         region_data.id(iter_region) = tmp_abbrv_to_id(tmp_region_abbrv);
%         region_data.Structure_name{iter_region} = registration_str.structure_table.name{registration_str.id_2_ind(region_data.id(iter_region))};
%     end
% end
% assert(all(isfinite(region_data.id)));
analysis_region_id_list = table_data.id;
analysis_region_name = table_data.name;
%% Settings for getting regional data
de_opt_str = struct;
de_opt_str.cube_stat_Q = true;
de_opt_str.node_feature_Q = false;
de_opt_str.link_feature_Q = false;
de_opt_str.vessel_graph_Q = false;
de_opt_str.depth_dependent_density_Q = false;
de_opt_str.merge_cc_Q = false;
de_opt_str.save_Q = false;
de_opt_str.save_fp = [];
merge_stack_name = 'all_stack';
allen_atlas_id = load('Allen_atlas_id.mat');
scatter_marker_type_list = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
num_maker_type = numel(scatter_marker_type_list);
if de_opt_str.merge_cc_Q
    num_cc = 1;
else
    num_cc = 2;
end

local_extrema_window_size_idx = 11; 
%% Collect all the region statistics in all the stack
region_data = cell(num_region, num_stack);
collect_data_tic = tic;
for iter_stack = 1 : num_stack
    for iter_region = 1 : num_region
        tmp_tic = tic;
        tmp_region_id = analysis_region_id_list(iter_region);
        tmp_region_name = strrep(analysis_region_name{iter_region}, ' ', '_');
        tmp_wb_data = wb_data_cell{iter_stack};
        tmp_region_data = fun_analysis_get_atlas_regional_data(tmp_wb_data, ...
            tmp_region_id, tmp_region_name, de_opt_str);
        %% Process region pO2 data
        tmp_pO2_data_field_name = setdiff(fieldnames(tmp_wb_data.cube_stat.pO2_data), ...
            {'dataset_name', 'stack', 'filepath', 'stack', 'skel_version', 'pO2_data_folder', ...
            'local_extrema_window_size'});
        if all(isfield(tmp_region_data, {'cc_1', 'cc_2'}))
            for iter_cc = 1 : tmp_region_data.num_cc
                tmp_cc_data = tmp_region_data.(sprintf('cc_%d', iter_cc));
                tmp_cc_data.local_cube_stat.num_cap = tmp_cc_data.local_cube_stat.link_cap_stat.num_data.length;
                tmp_cc_data.local_cube_stat.num_link = tmp_cc_data.local_cube_stat.link_all_stat.num_data.length;
                tmp_cc_data.local_cube_stat.node_density_mm3 = tmp_cc_data.local_cube_stat.node_stat.num_data.degree ./ (0.24 * linear_scaling_factor)^3;
                tmp_cc_data.local_cube_stat = rmfield(tmp_cc_data.local_cube_stat, {'mask_volume', 'num_surface_voxel'});
                
                for iter_field_name = 1 : numel(tmp_pO2_data_field_name)
                    tmp_field_name = tmp_pO2_data_field_name{iter_field_name};
                    tmp_cc_data.pO2_data.(tmp_field_name) = fun_structure_field_indexing(tmp_wb_data.cube_stat.pO2_data.(tmp_field_name), ...
                        tmp_cc_data.is_in_cc_cube_Q);
                end
                %% Compute local statistics
                tmp_cc_data.local_cube_stat.local_dt_mean = tmp_cc_data.pO2_data.local_dt_stat.mean;
                
                tmp_cc_data.local_cube_stat.local_pO2_mean = tmp_cc_data.pO2_data.local_pO2_stat.mean;
                
                tmp_cc_data.local_cube_stat.dt_lm_dt_v_mean = tmp_cc_data.pO2_data.dt_lm.dt_mean(:, local_extrema_window_size_idx);
                tmp_cc_data.local_cube_stat.dt_lm_pO2_v_mean = tmp_cc_data.pO2_data.dt_lm.pO2_mean(:, local_extrema_window_size_idx);
                tmp_cc_data.local_cube_stat.pO2_lm_dt_v_mean = tmp_cc_data.pO2_data.pO2_lm.dt_mean(:, local_extrema_window_size_idx);
                tmp_cc_data.local_cube_stat.pO2_lm_pO2_v_mean = tmp_cc_data.pO2_data.pO2_lm.pO2_mean(:, local_extrema_window_size_idx);
                
%                 tmp_cc_data.local_cube_stat.dt_lm_dt_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.dt_lm.v);
%                 tmp_cc_data.local_cube_stat.dt_lm_pO2_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.dt_lm.pO2_v);
%                 tmp_cc_data.local_cube_stat.pO2_lm_dt_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.pO2_lm.dt_v);
%                 tmp_cc_data.local_cube_stat.pO2_lm_pO2_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.pO2_lm.v);
                tmp_region_data.(sprintf('cc_%d', iter_cc)) = tmp_cc_data;
            end
        end
        region_data{iter_region, iter_stack} = tmp_region_data;
        toc(tmp_tic);
    end
end
fprintf('Finish collecting all the region data. Elapsed time is %f seconds.\n', toc(collect_data_tic));
%% Collect data
% 1. Capillary length density in 240-cube
% 2. Capillary segment length distribution PDF
% 3. Capillary-tissue maximum distance
min_in_cc_ratio = 0.75;
min_in_mask_ratio = 1;
min_cap2vsl_vol_ratio = 0.0;
fn_postfix = sprintf('iccr_%.2f_ibmr_%.2f_mc2vr_%.2f', min_in_cc_ratio, ...
    min_in_mask_ratio, min_cap2vsl_vol_ratio);
plot_data = struct;
[ plot_data.vessel_length_density, plot_data.capillary_length_density, ...
    plot_data.avg_dt_lm, plot_data.avg_d, ...
    plot_data.avg_u, plot_data.avg_u_lm, ...
    plot_data.cap_FA, plot_data.cap_FAz] = deal(cell(num_cc, num_stack, num_region));
for iter_region = 1 : num_region
    for iter_stack = 1 : num_stack
        tmp_region_data = region_data{iter_region, iter_stack};
        if all(isfield(tmp_region_data, {'cc_1', 'cc_2'}))
            for iter_cc = 1 : num_cc
                tmp_cc_data = tmp_region_data.(sprintf('cc_%d', iter_cc));
                tmp_240_cube_selected_Q = (tmp_cc_data.local_cube_stat.cube_in_brain_mask_ratio >= min_in_mask_ratio & ...
                    tmp_cc_data.local_cube_stat.in_cc_vol_f >= min_in_cc_ratio) & ...
                    tmp_cc_data.local_cube_stat.cap2vsl_vol_fraction >= min_cap2vsl_vol_ratio;
                plot_data.vessel_length_density{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.local_cube_stat.link_length_density_m_mm3(tmp_240_cube_selected_Q);
                plot_data.capillary_length_density{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.local_cube_stat.capillary_length_density_m_mm3(tmp_240_cube_selected_Q);
                
                plot_data.cap_FA{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.local_cube_stat.wb_ai_cap_vw.fractional_anisotropy(tmp_240_cube_selected_Q);
                plot_data.cap_FAz{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.local_cube_stat.wb_ai_cap_vw.fa_z(tmp_240_cube_selected_Q);            
                
                plot_data.avg_d{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.pO2_data.local_dt_stat.mean(tmp_240_cube_selected_Q);
                
                plot_data.avg_dt_lm{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.local_cube_stat.dt_lm_dt_v_mean(tmp_240_cube_selected_Q);
                
                plot_data.avg_u{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.pO2_data.local_pO2_stat.mean(tmp_240_cube_selected_Q);
                plot_data.avg_u_lm{iter_cc, iter_stack, iter_region} = ...
                    tmp_cc_data.local_cube_stat.pO2_lm_pO2_v_mean(tmp_240_cube_selected_Q);
            end
        end
    end
end
fprintf('Finish collecting data for visualization\n');
%% Compute the histogram
capillary_hist = struct;
compute_hist_field_name = fieldnames(plot_data);
percentile_list = [0, 0.1, 1, 2.5, 5, 10, 25, 50, 75, 90, 95, 97.5, 99, 99.9, 100];
%%
is_valid_Q = true(num_region, 1);
for iter_field = 1 : numel(compute_hist_field_name)
    tmp_field_name = compute_hist_field_name{iter_field};
    tmp_data = plot_data.(tmp_field_name);
%     tmp_prctile_v = cellfun(@(x) prctile(x, percentile_list), ...
%         tmp_data, 'UniformOutput', false);    
    
    tmp_mean = cellfun(@(x) mean(x, 'omitnan'), tmp_data);
%     tmp_median = cellfun(@median, tmp_data);
    tmp_std = cellfun(@(x) std(x, 'omitnan'), tmp_data);
    
    tmp_region_mean = mean(reshape(tmp_mean, 2 * num_stack, num_region), 1);
    tmp_region_std = mean(reshape(tmp_std, 2 * num_stack, num_region), 1);
    tmp_region_mean_std = std(reshape(tmp_mean, 2 * num_stack, num_region), 1, 1);
    tmp_region_std_std = std(reshape(tmp_std, 2 * num_stack, num_region), 1, 1);
    
    is_valid_Q = is_valid_Q & isfinite(tmp_region_mean');
    
    table_data.(sprintf('%s_mean_mean', tmp_field_name)) = tmp_region_mean';
    table_data.(sprintf('%s_mean_std', tmp_field_name)) = tmp_region_mean_std';
    table_data.(sprintf('%s_std_mean', tmp_field_name)) = tmp_region_std';
    table_data.(sprintf('%s_std_std', tmp_field_name)) = tmp_region_std_std';
end

num_cube_per_cc = cellfun(@numel, plot_data.vessel_length_density);
num_cube_per_cc = reshape(num_cube_per_cc, 2 * num_stack, num_region);
table_data.num_cube_mean = mean(num_cube_per_cc, 1, 'omitnan').';
table_data.num_cube_std = std(num_cube_per_cc, 1, 1, 'omitnan').';
%% 
output_filepath = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'all_stack'), ...
    'paper', sprintf('Allen_atlas_regional_data_%s.csv', fn_postfix));
output_folder = fileparts(output_filepath);
if ~isfolder(output_folder)
    mkdir(output_folder);
end
writetable(table_data, output_filepath);
output_filepath = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'all_stack'), ...
    'paper', sprintf('Allen_atlas_regional_data_%s_valid.csv', fn_postfix));

writetable(table_data(is_valid_Q, :), output_filepath);