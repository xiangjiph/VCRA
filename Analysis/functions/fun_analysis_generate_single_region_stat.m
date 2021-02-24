function exit_code = fun_analysis_generate_single_region_stat(dataset_name, stack_list, region_atlas_id)

exit_code = 1;
compute_region_stat_Q = false;
compute_240_cube_stat_Q = true;
% Extract field list
compute_240_cube_stat_field_list = {'mean'};
num_240_cube_stat_field = numel(compute_240_cube_stat_field_list);
cc_in_240_cube_data_field_list = {'link_all', 'link_cap'};
num_240_cube_stat_data = numel(cc_in_240_cube_data_field_list);

compute_network_properties_Q = true;
single_cc_Q = true;
single_stack_Q = true;
merged_stack_Q = true;

num_stack = numel(stack_list);
%% parameters
max_cap_r_um = 3.5;
min_in_mask_f = 0.95;
min_in_cc_f = 0.5;
%% Load the data from all the stack and hemispheres
region_data = fun_analysis_load_single_region_stat_data(dataset_name, stack_list, region_atlas_id);
region_name_abv = region_data.region_name_abv;
region_name = region_data.region_name;
valid_data_Q = region_data.valid_Q;
stack_data = region_data.data;
%% reorganize the stack data
str_info = struct;
str_info.dataset_name = dataset_name;
str_info.stack_list = stack_list;
str_info.cc_list = {1, 2};
str_info.num_cc = 2;
str_info.region_name = region_name;
str_info.region_name_abv = region_name_abv;
str_info.dim_name = {'cc', 'stack'};
str_info.dim_name_list = {{'cc_1', 'cc_2'}, stack_list};

if all(valid_data_Q)
    if compute_region_stat_Q
        region_stat = struct;
        region_stat.info = str_info;
        region_stat.info.str_name = 'region_stat';
        [region_stat.link_all, region_stat.link_cap, region_stat.node] = deal(cell(2, num_stack));
    end
    if compute_240_cube_stat_Q
        cube_feature_stat = struct;
        cube_feature_stat.info = str_info;
        cube_feature_stat.info.str_name = 'region_240_cubes_feature_stat_mean';
        [cube_feature_stat.link_all, cube_feature_stat.link_cap, cube_feature_stat.node] = deal(cell(2, num_stack));
    end
    if compute_network_properties_Q
        cube_stat = struct;
        cube_stat.info = str_info;
        cube_stat.info.str_name = 'region_240_cubes_stat';
        cube_stat.network_prop = deal(cell(2, num_stack));
    end
    for iter_stack = 1 : num_stack
        for iter_cc = 1 : 2
            tmp_cc_data = stack_data{iter_stack}.(sprintf('cc_%d', iter_cc));
            if compute_region_stat_Q
                region_stat.link_all{iter_cc, iter_stack} = tmp_cc_data.link_features;
                region_stat.link_cap{iter_cc, iter_stack} = tmp_cc_data.link_features(tmp_cc_data.link_features.dt_median <= max_cap_r_um, :);
                region_stat.node{iter_cc, iter_stack} = tmp_cc_data.node_features;
            end
            
            if compute_network_properties_Q
                tmp_network_data = rmfield(tmp_cc_data.local_cube_stat, {'link_all', 'link_cap', ...
                    'node', 'ai_all_vw_vec', 'ai_cap_vw_vec', 'ai_noncap_vw_vec', 'grid_sub'});
            end
            % Preselect internal cubes
            tmp_selected_Q = true(size(tmp_network_data.volume_ratio));
            if isfield(tmp_network_data, 'in_mask_vol_f')
                tmp_selected_Q = tmp_selected_Q & tmp_network_data.in_mask_vol_f > min_in_mask_f;
            end
            if isfield(tmp_network_data, 'in_cc_vol_f')
                tmp_selected_Q = tmp_selected_Q & tmp_network_data.in_cc_vol_f > min_in_cc_f;
            end
            
            if compute_240_cube_stat_Q
                for iter_data = 1 : num_240_cube_stat_data
                    tmp_stat_data_name = cc_in_240_cube_data_field_list{iter_data};
                    tmp_stat_str = struct;
                    for iter_stat_field = 1 : num_240_cube_stat_field
                        tmp_stat_field_name = compute_240_cube_stat_field_list{iter_stat_field};
                        tmp_cc_feature_data = tmp_cc_data.local_cube_stat.(tmp_stat_data_name);
                        
                        tmp_cc_feature_data = tmp_cc_feature_data(tmp_selected_Q);
           
                        tmp_stat_str.(tmp_stat_field_name) = ...
                            fun_analysis_get_fields_in_cell_array_str(tmp_cc_feature_data, ...
                            compute_240_cube_stat_field_list{iter_stat_field});
                    end
                    tmp_stat_str = tmp_stat_str.mean;
                    if ~isempty(tmp_stat_str)
                        tmp_stat_str = struct2table(tmp_stat_str); % hard code here
                        cube_feature_stat.(tmp_stat_data_name){iter_cc, iter_stack} = ...
                            tmp_stat_str;
                    else
                        fprintf('The 240 cube statistics is empty. Return\n');
                        return;
                    end
                end
            end            
            if compute_network_properties_Q
                % Convert to table here. isfield DOES NOT WORK for table!
                tmp_network_data = struct2table(tmp_network_data);
                tmp_network_data = tmp_network_data(tmp_selected_Q, :);
                cube_stat.network_prop{iter_cc, iter_stack} = tmp_network_data;
            end            
        end
    end
else
    fprintf('Not all the stacks are valid\n');
    exit_code = 1;
    return;
end
%% Define plot settings
if compute_region_stat_Q
    link_plot_setting = {'length', 'dt_median', 'surface_area', 'volume', ...
        'nearest_tissue_radius', 'capillary_branching_order', 'shortest_loop_length', ...
        'shortest_loop_geodesic_length', 'dist_to_nearest_noncapillary_mean', ...
        'straightness'};    
    region_stat.plot_setting.link_all = link_plot_setting;
    region_stat.plot_setting.link_cap = link_plot_setting;
    
    node_plot_setting = {'nearest_node_dist', 'path_to_nearest_neighbor_geodesic_length', ...
        'degree', 'link_length_min', 'link_length_max', 'link_length_median'};
    region_stat.plot_setting.node = node_plot_setting;
    % Plot histogram for each cc
    if single_cc_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(region_stat, []);
    end
    if single_stack_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(region_stat, 1);
    end
    if merged_stack_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(region_stat, [1, 2]);
    end
end
%% Local 240 cube cc feature properties
if compute_240_cube_stat_Q
    link_cc_stat_setting = {'length', 'nearest_tissue_dt_mean', 'nearest_tissue_dt_median', ...
        'nearest_tissue_dt_max', 'nearest_tissue_radius'};
    cube_feature_stat.plot_setting.link_all = link_cc_stat_setting;
    cube_feature_stat.plot_setting.link_cap = link_cc_stat_setting;

    if single_cc_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(cube_feature_stat, []);
    end
    if single_stack_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(cube_feature_stat, 1);
    end
    if merged_stack_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(cube_feature_stat, [1, 2]);
    end    
end
%% Local 240 network properties
if compute_network_properties_Q
    network_plot_setting = {'surface_area_mm2mm3', 'volume_ratio', 'length_density_m_mm3', ...
        'length_density_cap_m_mm3', ...
        'ai_all_vw_fa', 'ai_all_vw_fa_z', 'ai_all_vw_min2max_z', 'ai_all_vw_svrmax', 'ai_all_vw_svd1_z', ...
        'ai_cap_vw_fa', 'ai_cap_vw_fa_z', 'ai_cap_vw_svrmax', 'ai_cap_vw_min2max_z', 'ai_cap_vw_svd1_z'};
    
    cube_stat.plot_setting.network_prop = network_plot_setting;
    if single_cc_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(cube_stat, []);
    end
    if single_stack_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(cube_stat, 1);
    end
    if merged_stack_Q
        exit_code = fun_analysis_compute_histogram_of_region_stat(cube_stat, [1,2]);
    end
end
end
