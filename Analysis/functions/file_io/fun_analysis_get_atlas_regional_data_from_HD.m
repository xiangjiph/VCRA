function region_stat = fun_analysis_get_atlas_regional_data_from_HD(registration_str, grid_info, wb_mask_str, ...
    skeleton_version, mask_version, id_to_merge_subregion, merged_str_name, save_cc_featureQ)
if nargin < 6
    merged_str_name = [];
    save_cc_featureQ = false;
elseif nargin < 7
    save_cc_featureQ = false;
end

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end
max_cap_r_um = 3.5;
%% Parse input
dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
image_grid_version = grid_info.version;

wb_mask = wb_mask_str.dt > 0;
wb_mask_dt = wb_mask_str.dt;
wb_mask_ds_ratio = wb_mask_str.ds_ratio;

if ~isempty(merged_str_name)
    name_to_merge_subregion = merged_str_name;
else
    name_to_merge_subregion = registration_str.label_2_name(id_to_merge_subregion);
end
name_to_save = strrep(name_to_merge_subregion, ' ', '_');
%% Save information
region_stat = struct;
region_stat.dataset_name = dataset_name;
region_stat.stack = stack;
region_stat.registration_data_fp = registration_str.info.fp;
region_stat.grid_version = image_grid_version;
region_stat.feature_version = mask_version;
region_stat.skeleton_version = skeleton_version;
region_stat.structure_name = name_to_merge_subregion;
region_stat.folder_path = fullfile(DataManager.fp_metadata_folder(dataset_name, stack), ...
    'region_data');
region_stat.filepath = fullfile(region_stat.folder_path, ...
    sprintf('%s_%s_region_data_%s.mat', dataset_name, stack, name_to_save));
%% Merge atals region by annotation path depth
[subregion_cc, subregion_mask] = fun_registration_get_region_cc(registration_str, id_to_merge_subregion);

for cc_idx = 1 : subregion_cc.NumObjects
    tmp_tic = tic;
    cc_data = struct;
    cc_data.structure_name = name_to_merge_subregion;
    cc_data.array_size = registration_str.array_size;
    cc_data.cc_ind = subregion_cc.PixelIdxList{cc_idx};
    cc_data.voxel_size = registration_str.info.voxel_size_um;
    fprintf('Processing %s connected component %d/%d.\n', name_to_merge_subregion, cc_idx, subregion_cc.NumObjects);
    cc_mask = (subregion_mask == cc_idx);
    %% Determine the 240-cubes that cover the structure
    in_cc_volume_fraction = fun_analysis_get_bbox_in_mask_vol_ratio(cc_mask, grid_info.bbox_xyz_mmxx_list ./ ...
        [registration_str.info.voxel_size_um, registration_str.info.voxel_size_um]);
    inside_cube_Q = in_cc_volume_fraction > 0;
    inside_cube_grid_sub = grid_info.bbox_grid_sub_list(inside_cube_Q, :);
    if ~any(inside_cube_Q)
        break;
    end
    %% Load all the skeleton inside the selected region
    fprintf('Loading region skeleton data to reconstruct the graph\n');
    region_graph = fun_graph_construct_graph_in_region_by_grid_sub_from_HD(dataset_name, stack, skeleton_version, inside_cube_grid_sub);
    fprintf('Loading region links and nodes features\n');
    [region_feature, local_cube_stat] = fun_analysis_load_features_in_region_by_grid_sub_from_HD(dataset_name, stack, mask_version, inside_cube_grid_sub);
    
    local_cube_stat.in_mask_vol_f = fun_analysis_get_bbox_in_mask_vol_ratio(wb_mask, grid_info.bbox_xyz_mmxx_list(inside_cube_Q, :) ./ wb_mask_ds_ratio);
    local_cube_stat.in_cc_vol_f = in_cc_volume_fraction(inside_cube_Q);

    fprintf('Get the link and node features from the loaded tabel\n');
    region_graph.link.features = fun_analysis_get_link_features_from_local_stat_table(region_graph, region_feature.link);
    region_graph.node.features = fun_analysis_get_node_features_from_local_stat_table(region_graph, region_feature.node);
    
    cc_data.local_cube_stat = local_cube_stat;
    %% Look at the regional statistics - link
    % Local reconstruction mask of the connected component - for
    % determining the fraction of links inside the structure mask
    
    % Find the link that are in the cc local mask
    fprintf('Compute the fraction of link inside the structure\n');
    link_in_cc_ratio = nan(region_graph.link.num_cc, 1);
    link_ind_ds = cell(region_graph.link.num_cc, 1);
    for iter_link = 1 : region_graph.link.num_cc
        tmp_ind = region_graph.link.cc_ind{iter_link};
        % Convert the coordiante to the local coordinate of the cc_local_mask
        tmp_sub = fun_ind2sub(region_graph.num.mask_size, tmp_ind);
        tmp_sub = round(tmp_sub ./ registration_str.info.voxel_size_um);
        tmp_sub = min(max(1, tmp_sub), registration_str.array_size);
        tmp_ind_rs = sub2ind(registration_str.array_size, tmp_sub(:, 1), tmp_sub(:, 2), ...
            tmp_sub(:, 3));
        link_in_cc_ratio(iter_link) = nnz(cc_mask(tmp_ind_rs)) / numel(tmp_ind_rs);
        link_ind_ds{iter_link} = tmp_ind_rs;
    end
    % Fraction of links that are in the connected component but not valid in
    % the transferred link features
    link_has_featureQ = ~isnan(region_graph.link.features.length);
    num_link_wo_feature_in_cc = nnz(~link_has_featureQ & link_in_cc_ratio > 0);
    fprintf('\t%d/%d links in the connected component cannot find features from the local analysis table\n', num_link_wo_feature_in_cc, numel(link_has_featureQ));
    %% Compute depth-dependence on the capillary density
    str_cc_dt = wb_mask_dt(cc_mask);
    cc_volume = numel(str_cc_dt) * prod(registration_str.info.voxel_size_um);
    cc_data.volume = cc_volume;
    
    fprintf('%f of the structure mask is inside the brain mask\n', nnz(str_cc_dt)/(numel(str_cc_dt)));
    % Select the part of the structure mask that are inside the
    % brain mask
    str_cc_dt = str_cc_dt(str_cc_dt > 0);
    if isempty(str_cc_dt)
        break;
    end
    bin_dt_edge = prctile(str_cc_dt(:), 0.5) : 50 : prctile(str_cc_dt(:), 99.5);
    str_cc_dt_bin = histcounts(str_cc_dt, bin_dt_edge);
    str_cc_dt_bin = str_cc_dt_bin .* (wb_mask_ds_ratio ^ 3); % Voxel size is (16 um)^3
    cc_data.dt.bin_edge = bin_dt_edge;
    cc_data.dt.cc_voxel_count = str_cc_dt_bin;
    
    is_in_cc_Q = link_in_cc_ratio > 0;
    is_capillary_link_in_cc_Q = (region_graph.link.features.dt_median <= max_cap_r_um & is_in_cc_Q);
    
    region_ves_dt = wb_mask_dt(cat(1, link_ind_ds{is_in_cc_Q}));
    region_cap_dt = wb_mask_dt(cat(1, link_ind_ds{is_capillary_link_in_cc_Q}));
    fprintf('%f of the skeleton voxel is inside the brain mask\n', nnz(region_cap_dt)/(numel(region_cap_dt)));
    region_cap_dt = region_cap_dt(region_cap_dt > 0);
    region_ves_dt = region_ves_dt(region_ves_dt > 0);
    
    region_cap_dt_bin = histcounts(region_cap_dt, bin_dt_edge);
    region_ves_dt_bin = histcounts(region_ves_dt, bin_dt_edge);
    
    cc_data.dt.cap_voxel_count = region_cap_dt_bin;
    cc_data.dt.ves_voxel_count = region_ves_dt_bin;
    % Compute the scaling factor to convert number of skeleton to
    % length
    ratio_numskel_to_length = region_graph.link.features.length ./ region_graph.link.num_voxel_per_cc;
    mean_ratio = mean(ratio_numskel_to_length, 'omitnan');
    cc_data.dt.numskel2length_ratio = mean_ratio;
    
    cc_data.dt.cap_length_density = region_cap_dt_bin .* mean_ratio ./ (1e6 * str_cc_dt_bin / 1e9);
    cc_data.dt.ves_length_density = region_ves_dt_bin .* mean_ratio ./ (1e6 * str_cc_dt_bin / 1e9);
    cc_data.dt.bin_val = movmean(bin_dt_edge, 2, 'Endpoints', 'discard');
    %% Regional statistics
    cc_link_features = region_graph.link.features(is_in_cc_Q, :);
    
    cc_link_features.in_cc_ratio = link_in_cc_ratio(is_in_cc_Q);
    cc_link_features.length_in_cc = cc_link_features.length .* cc_link_features.in_cc_ratio;
    cc_link_stat = fun_analysis_get_graph_feature_basic_stat(cc_link_features);
    if save_cc_featureQ
        cc_data.link_features = cc_link_features;
    end
    cc_data.link_stat = cc_link_stat;
    %% Regional statistics - node
    node_ind_ds = fun_analysis_get_scaled_cc_ind({region_graph.node.features.global_ind}, ...
        region_graph.num.mask_size, 1 ./ registration_str.info.voxel_size_um, registration_str.array_size);
    node_ind_ds = node_ind_ds{:};
    node_in_cc_Q = cc_mask(node_ind_ds);
    cc_node_features = region_graph.node.features(node_in_cc_Q, :);
    cc_node_stat = fun_analysis_get_graph_feature_basic_stat(cc_node_features);
    if save_cc_featureQ
        cc_data.node_features = cc_node_features;
    end
    cc_data.node_stat = cc_node_stat;
    %% Record the cc information
    region_stat.(sprintf('cc_%d', cc_idx)) = cc_data;
    fprintf('Finish processing %s connected component %d. Elapsed time is %f seconds\n',...
        name_to_merge_subregion, cc_idx,  toc(tmp_tic));
end
%% Save data
if ~isfolder(region_stat.folder_path)
    mkdir(region_stat.folder_path);
end
save(region_stat.filepath, '-struct', 'region_stat', '-v7.3');
end