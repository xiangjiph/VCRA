function exit_code = fun_vis_generate_single_region_stat_visualization(dataset_name, stack_list, region_atlas_id)

persistent registration_str DataManager;
if isempty(registration_str)
    registration_str = load('Allen_atlas.mat');
end
if isempty(DataManager)
    DataManager = FileManager;
end

compute_region_stat_Q = true;
compute_network_properties_Q = true;
plot_cc_Q = false;
plot_stack_Q = true;
plot_merged_stack_Q = true;
%% parameters
max_cap_r_um = 3.5;
%% Load the data from all the stack and hemispheres
if isnumeric(region_atlas_id)
    region_list_ind = full(registration_str.id_2_ind(region_atlas_id));
    if region_list_ind > 0
        region_name = registration_str.structure_table.name(region_list_ind);
        region_name = region_name{1};
        region_name_abv = registration_str.structure_table.acronym(region_list_ind);
        region_name_abv = region_name_abv{1};
    else
       error('The structure atlas id does not exist'); 
    end
else
    region_name = region_atlas_id;
    region_name_abv = region_atlas_id;    
end
num_stack = numel(stack_list);
stack_data = cell(num_stack, 1);
valid_data_Q = false(num_stack, 1);
for iter_stack = 1 : num_stack
    tmp_stack = stack_list{iter_stack};
    region_data_folder = fullfile(DataManager.fp_metadata_folder(dataset_name, ...
        tmp_stack), 'region_data');
    region_data_filepath = fullfile(region_data_folder, ...
        sprintf('%s_%s_region_data_%s.mat', dataset_name, tmp_stack, ...
        strrep(region_name, ' ', '_')));
    if isfile(region_data_filepath)
        tmp_tic = tic;
        stack_data{iter_stack} = load(region_data_filepath);
        fprintf('Finish loading %s. Elapsed time is %f seconds\n', ...
            region_data_filepath, toc(tmp_tic));
        if isfield(stack_data{iter_stack}, 'cc_1') && isfield(stack_data{iter_stack}, 'cc_2')
            valid_data_Q(iter_stack) = true;
        end
    else
        fprintf('File %s does not exist\n', region_data_filepath);
    end
end
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
    %     cube_feature_stat = struct;
    %     cube_feature_stat.info = str_info;
    %     cube_feature_stat.info.str_name = 'region_240_cubes_feature_stat';
    %     [cube_feature_stat.link_all, cube_feature_stat.link_cap, cube_feature_stat.node] = deal(cell(2, num_stack));
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
            %             cube_feature_stat.link_all{iter_cc, iter_stack} = tmp_cc_data.local_cube_stat.link_all;
            %             cube_feature_stat.link_cap{iter_cc, iter_stack} = tmp_cc_data.local_cube_stat.link_cap;
            %             cube_feature_stat.node{iter_cc, iter_stack} = tmp_cc_data.local_cube_stat.node;
            if compute_network_properties_Q
                cube_stat.network_prop{iter_cc, iter_stack} = rmfield(tmp_cc_data.local_cube_stat, {'link_all', 'link_cap', ...
                    'node', 'ai_all_vw_vec', 'ai_cap_vw_vec', 'ai_noncap_vw_vec', 'grid_sub'});
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
    link_plot_setting = cell(4, 0);
    % Format: variable name, label, xscale, yscale
    link_plot_setting(:, end+1) = {'length', 'Vessel segment length(\mum)', 'linear', 'log'};
    link_plot_setting(:, end+1) = {'dt_median', 'Vessel radius(\mum)', 'linear', 'log'};
    link_plot_setting(:, end+1) = {'surface_area', 'Vessel segment surface area(\mum^2)', 'linear', 'log'};
    link_plot_setting(:, end+1) = {'volume', 'Vessel segment volume(\mum^3)', 'linear', 'log'};
    
    link_plot_setting(:, end+1) = {'nearest_tissue_radius', 'Nearest tissue radius(\mum)', 'linear', 'log'};
    
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
    network_plot_setting(:, end+1) = {'surface_area_mm2mm3', 'Vessel surface area density(mm^2/mm^3)', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'volume_ratio', 'Vessel volume ratio', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'length_density_m_mm3', 'Vessel length density(m/mm^3)', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'length_density_cap_m_mm3', 'Capillary length density(m/mm^3)', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_all_vw_fa', 'Vessel fractional anisotropy', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_all_vw_fa_z', 'Vessel fractional anisotropy score', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_all_vw_min2max_z', 'Vessel isotropy score', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_all_vw_svrmax', 'Vessel orientation max PCA ratio', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_all_vw_svd1_z', 'Vessel anisotropy score', 'linear', 'linear'};
    
    network_plot_setting(:, end+1) = {'ai_cap_vw_fa', 'Capillary fractional anisotropy', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_cap_vw_fa_z', 'Capillary fractional anisotropy score', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_cap_vw_svrmax', 'Capillary orientation max PCA ratio', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_cap_vw_min2max_z', 'Capillary isotropy score', 'linear', 'linear'};
    network_plot_setting(:, end+1) = {'ai_cap_vw_svd1_z', 'Capillary anisotropy score', 'linear', 'linear'};


    % Preselect internal cubes
    min_in_mask_f = 0.95;
    min_in_cc_f = 0.5;
    for iter_cell = 1 : numel(cube_stat.network_prop)
        tmp_data = cube_stat.network_prop{iter_cell};
        tmp_selected_Q = true(size(tmp_data.volume_ratio));
        if isfield(tmp_data, 'in_mask_vol_f')
            tmp_selected_Q = tmp_selected_Q & tmp_data.in_mask_vol_f > min_in_mask_f;
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
