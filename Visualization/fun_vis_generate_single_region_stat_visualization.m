function exit_code = fun_vis_generate_single_region_stat_visualization(stack_data)

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end

compute_region_feature_stat_Q = true;

compute_network_properties_Q = true;

compute_region_cube_feature_stat_Q = true;
% Preselect internal cubes
min_in_mask_f = 0.95;
min_in_cc_f = 0.5;

plot_cc_Q = false;
plot_stack_Q = true;
plot_merged_stack_Q = true;

extract_ai_name_list = {'wb_ai_all_vw', 'wb_ai_cap_vw'};
extract_ai_field_name_list = {'fractional_anisotropy', 'fa_z', 'fa_p', 'svd_value_ratio'};

%% Parse input
region_name = strrep(stack_data{1}.structure_name, '_', ' ');
region_name_abv = region_name; % Add abbreviation later
num_stack = numel(stack_data);
valid_data_Q = cellfun(@(x) isfield(x, 'cc_1') && isfield(x, 'cc_2'), stack_data);
%% reorganize the stack data
str_info = struct;
str_info.dataset_name = stack_data{1}.dataset_name;
str_info.stack_list = cellfun(@(x) x.stack, stack_data, 'UniformOutput', false);
str_info.cc_list = {1, 2};
str_info.num_cc = numel(str_info.cc_list);
str_info.region_name = region_name;
str_info.region_name_abv = region_name_abv;
str_info.dim_name = {'cc', 'stack'};
str_info.dim_name_list = {{'cc_1', 'cc_2'}, str_info.stack_list};

if all(valid_data_Q)
    if compute_region_feature_stat_Q
        region_stat = struct;
        region_stat.info = str_info;
        region_stat.info.str_name = 'region_stat';
        [region_stat.link_all, region_stat.link_cap, region_stat.node] = deal(cell(str_info.num_cc, num_stack));
    end
    if compute_region_cube_feature_stat_Q
        cube_feature_stat = struct;
        cube_feature_stat.info = str_info;
        cube_feature_stat.info.str_name = 'region_240_cubes_feature_stat';
        [cube_feature_stat.link_all, cube_feature_stat.link_cap, cube_feature_stat.node] = deal(cell(str_info.num_cc, num_stack));
    end
    if compute_network_properties_Q
        cube_stat = struct;
        cube_stat.info = str_info;
        cube_stat.info.str_name = 'region_240_cubes_stat';
        cube_stat.network_prop = deal(cell(str_info.num_cc, num_stack));
    end
    for iter_stack = 1 : num_stack
        for iter_cc = 1 : 2
            tmp_cc_data = stack_data{iter_stack}.(sprintf('cc_%d', iter_cc));
            if compute_region_feature_stat_Q
                region_stat.link_all{iter_cc, iter_stack} = tmp_cc_data.link_features;
                region_stat.link_cap{iter_cc, iter_stack} = tmp_cc_data.link_features(tmp_cc_data.link_features.is_large_vessel_Q == 0, :);
                region_stat.node{iter_cc, iter_stack} = tmp_cc_data.node_features;
            end
            if compute_region_cube_feature_stat_Q            
                cube_feature_stat.link_all{iter_cc, iter_stack} = tmp_cc_data.local_cube_stat.link_all_stat;
                cube_feature_stat.link_cap{iter_cc, iter_stack} = tmp_cc_data.local_cube_stat.link_cap_stat;
                cube_feature_stat.node{iter_cc, iter_stack} = tmp_cc_data.local_cube_stat.node_stat;
            end
            if compute_network_properties_Q
                tmp_network_prop_str = rmfield(tmp_cc_data.local_cube_stat, {'link_all_stat', 'link_cap_stat', ...
                    'node_stat', 'dataset_name', 'grid_version', 'mask_version', 'stack', ...
                    'wb_ai_all_lw', 'wb_ai_cap_lw'});
                % Extract fields                
                for iter_ai = 1 : numel(extract_ai_name_list)
                    tmp_ai_name = extract_ai_name_list{iter_ai};
                    for iter_ai_field = 1 : numel(extract_ai_field_name_list)
                        tmp_ai_field_name = extract_ai_field_name_list{iter_ai_field};
                        tmp_network_prop_str.(sprintf('%s_%s', tmp_ai_name, tmp_ai_field_name)) = ...
                            tmp_network_prop_str.(tmp_ai_name).(tmp_ai_field_name)(:, 1);
                    end
                end
                cube_stat.network_prop{iter_cc, iter_stack} = rmfield(tmp_network_prop_str, {'wb_ai_all_vw', 'wb_ai_cap_vw'});
            end
        end
    end
else
    fprintf('Not all the stacks are valid\n');
    exit_code = 1;
    return;
end
%% Define plot settings
if compute_region_feature_stat_Q
    link_plot_setting = cell(4, 0);
    % Format: variable name, label, xscale, yscale
    link_plot_setting(:, end+1) = {'length', 'Vessel segment length(\mum)', 'linear', 'log'};
    link_plot_setting(:, end+1) = {'dt_median', 'Vessel radius(\mum)', 'linear', 'log'};
    link_plot_setting(:, end+1) = {'surface_area', 'Vessel segment surface area(\mum^2)', 'linear', 'log'};
    link_plot_setting(:, end+1) = {'volume', 'Vessel segment volume(\mum^3)', 'linear', 'log'};
    
    link_plot_setting(:, end+1) = {'nearest_tissue_radius', 'Nearest tissue radius(\mum)', 'linear', 'log'};
    link_plot_setting(:, end+1) = {'nearest_tissue_dt_max', 'Nearest tissue maximum distance (\mum)', 'linear', 'linear'};
    
    link_plot_setting(:, end+1) = {'capillary_branching_order', 'Capillary branching order', 'linear', 'linear'};
    link_plot_setting(:, end+1) = {'shortest_loop_length', 'Length of shortest loop(\mum)', 'linear', 'linear'};
    link_plot_setting(:, end+1) = {'shortest_loop_geodesic_length', 'Number of edges', 'linear', 'linear'};
    
    link_plot_setting(:, end+1) = {'dist_to_nearest_noncapillary_mean', 'Distance to the nearest noncapillary(\mum)', 'linear', 'linear'};
    link_plot_setting(:, end+1) = {'straightness', 'Straightness', 'linear', 'linear'};
    
    region_stat.plot_setting.link_all = link_plot_setting;
    region_stat.plot_setting.link_cap = link_plot_setting;
    
    node_plot_setting = cell(4, 0);
    node_plot_setting(:, end+1) = {'nearest_node_dist', 'Distance to the nearest node(\mum)', 'linear', 'linear'};
    node_plot_setting(:, end+1) = {'path_to_nearest_neighbor_geodesic_length', 'Number of edges to the nearest node', 'linear', 'linear'};
    node_plot_setting(:, end+1) = {'degree', 'degree', 'linear', 'linear'};
    node_plot_setting(:, end+1) = {'link_length_min', 'Length of the shortest link(\mum)', 'linear', 'linear'};
    node_plot_setting(:, end+1) = {'link_length_max', 'Length of the longest link(\mum)', 'linear', 'linear'};
    node_plot_setting(:, end+1) = {'link_length_median', 'Median length of the links(\mum)', 'linear', 'linear'};
    region_stat.plot_setting.node = node_plot_setting;
    % Plot histogram for each cc
    if plot_cc_Q
        exit_code = fun_vis_histogram_of_region_stat(region_stat, []);
    end
    if plot_stack_Q
        exit_code = fun_vis_histogram_of_region_stat(region_stat, 1);
    end
    if plot_merged_stack_Q
        exit_code = fun_vis_histogram_of_region_stat(region_stat, [1, 2]);
    end
end
%% Local 240 network properties
if compute_network_properties_Q
    network_plot_setting = cell(4, 0);
    
    network_plot_setting(:, end+1) = {'link_surface_area_density_mm2_mm3', 'Vessel surface area density(mm^2/mm^3)', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'link_volume_density', 'Vessel volume ratio', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'link_length_density_m_mm3', 'Vessel length density(m/mm^3)', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'capillary_length_density_m_mm3', 'Capillary length density(m/mm^3)', 'linear', 'linear'};
    
    network_plot_setting(:, end+1) = {'wb_ai_all_vw_fractional_anisotropy', 'Vessel fractional anisotropy', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'wb_ai_all_vw_fa_z', 'Vessel FA test z-score', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'wb_ai_all_vw_fa_p', 'Vessel FA test p-value', 'log', 'linear'};
    network_plot_setting(:, end+1) = {'wb_ai_all_vw_svd_value_ratio', 'Vessel PCV', 'linear', 'linear'};
    
    network_plot_setting(:, end+1) = {'wb_ai_cap_vw_fractional_anisotropy', 'Capillary fractional anisotropy', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'wb_ai_cap_vw_fa_z', 'Capillary FA test z-score', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'wb_ai_cap_vw_fa_p', 'Capillary FA test p-value', 'log', 'linear'};
    network_plot_setting(:, end+1) = {'wb_ai_cap_vw_svd_value_ratio', 'Capillary PCV', 'linear', 'linear'};

    network_plot_setting(:, end+1) = {'cap2vsl_length_fraction', 'Capillary-vessel length density ratio', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'cap2vsl_vol_fraction', 'Capillary-vessel volume density ratio', 'linear', 'linear'};

    for iter_cell = 1 : numel(cube_stat.network_prop)
        tmp_data = cube_stat.network_prop{iter_cell};
        tmp_selected_Q = true(size(tmp_data.cube_in_brain_mask_ratio));
        if isfield(tmp_data, 'in_mask_vol_f')
            tmp_selected_Q = tmp_selected_Q & tmp_data.cube_in_brain_mask_ratio > min_in_mask_f;
        end
        if isfield(tmp_data, 'in_cc_vol_f')
            tmp_selected_Q = tmp_selected_Q & tmp_data.in_cc_vol_f > min_in_cc_f;
        end
        % Convert to table here. isfield DOES NOT WORK for table!
        tmp_data = struct2table(tmp_data);
        tmp_data = tmp_data(tmp_selected_Q, :);
        cube_stat.network_prop{iter_cell} = tmp_data;
    end
    cube_stat.plot_setting.network_prop = network_plot_setting;
    if plot_cc_Q
        exit_code = fun_vis_histogram_of_region_stat(cube_stat, []);
    end
    if plot_stack_Q
        exit_code = fun_vis_histogram_of_region_stat(cube_stat, 1);
    end
    if plot_merged_stack_Q
        exit_code = fun_vis_histogram_of_region_stat(cube_stat, [1,2]);
    end
end
end
