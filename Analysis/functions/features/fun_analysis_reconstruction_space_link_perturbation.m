function dtf = fun_analysis_reconstruction_space_link_perturbation(mask_labeled, vessel_graph, ...
    connectivity_graph, link_label, downsample_scale, nonmax_win_size_um)

if nargin < 4
    link_label = vessel_graph.link.label;
    downsample_scale = 0.5;
    nonmax_win_size_um = 40; % Size of the nonmaximum suppression window. The actual phsyics size need to be rescaled. 
elseif nargin < 5
    downsample_scale = 0.5;
    nonmax_win_size_um = 40;
elseif nargin < 6
    nonmax_win_size_um = 40;
end
num_label = numel(link_label);


tissue_dt_prctile = [0, 0.1, 1, 2.5, 5 : 5 : 95, 97.5, 99, 99.9, 100];
num_prctile = numel(tissue_dt_prctile);
%% Initialization
% Unperturbed tissue space distance transform properties
dtf = struct;
[dtf.nearest_tissue_dt_mean, dtf.nearest_tissue_dt_max, ...
    dtf.nearest_tissue_dt_median, dtf.nearest_tissue_volume, ...
    dtf.recon_mask_num_voxel] = deal(nan(num_label, 1));
dtf.tissue_dt_prctile = nan(num_prctile, num_label);
% Perturbed tissue space distance transform properties
[dtf.nb_lk_vol_r, dtf.nb_lk_bch_od, dtf.nb_lk_label] = deal(cell(num_label, 1));
[dtf.num_nb_lk, ...
    dtf.nb_lk_vol_r_max, ...
    dtf.nb_lk_bch_od_max, ...
    dtf.nb_lk_bch_od_max_vol_r, dtf.nb_lk_vol_r_max_bch_od, ...
    dtf.nb_lk_bch_od_1_vol_r, ...
    dtf.dist_rm_lk_2_nlm, dtf.nlm_v_af_rm, ...
    dtf.lk_dt_af_rm_mean, dtf.lk_dt_af_rm_median, dtf.lk_dt_af_rm_max, ...
    dtf.pt_vol_dt_max, dtf.dist_rm_lk_2_pt_max_dt, ...
    dtf.pt_vol_dt_mean, dtf.pt_vol_dt_median, ...
    dtf.up_ts_2_p_ts_dt_mean, ...
    dtf.diff_dp_d_mean] = deal(nan(num_label, 1));
if num_label == 0
    return;
end
%% Parameters
bbox_expand_length = 3;
max_valid_node_label = connectivity_graph.numnodes;
% Determine node in the self-loop, which have been excluded in the
% connectivity graph
% self_loop_link_label = find(vessel_graph.link.connected_node_label(:, 1) == ...
%     vessel_graph.link.connected_node_label(:, 2) & vessel_graph.link.connected_node_label(:, 1) > 0);
% self_loop_node_label = vessel_graph.link.connected_node_label(self_loop_link_label, :);
%% Compute distance transform properties before perturbation: 
if downsample_scale ~= 1
    % Only nearest neighbor interpolation is applicable, since we don't
    % want to change the value of the labels
    % Downsampling the mask would dramatically save the memory and the
    % computation time. In the test block, the distance are overestimated
    % by about 3% and the volume is underestimated.
    mask_labeled = imresize3(mask_labeled, downsample_scale, 'Method', 'nearest');
end
nonmax_win_size = nonmax_win_size_um * downsample_scale;

image_size = size(mask_labeled);
im_dims = nnz(image_size ~= 1);
num_vol_vxl = numel(mask_labeled);
% Compute distance between the tissue to the nearest vessel
[tissue_space_dt, tissue_nearest_mask_ind] = bwdist(mask_labeled~=0);
tissue_space_dt = tissue_space_dt ./ downsample_scale;
valid_tissue_voxel_ind = find(tissue_space_dt > 0);
if num_vol_vxl  < 65535
    valid_tissue_voxel_ind = uint16(valid_tissue_voxel_ind);
elseif num_vol_vxl  < 4294967295
    valid_tissue_voxel_ind = uint32(valid_tissue_voxel_ind);
end
% The mask voxel indices (before cropping) that are closest to the space point
% valid_nearest_mask_ind_list = tissue_nearest_mask_ind(valid_tissue_voxel_ind);
valid_tissue_space_dt_list = tissue_space_dt(valid_tissue_voxel_ind);
clearvars tissue_space_dt
% The nearest label is from the graph that reconstruct the labeled matrix,
% not the graph that use for internal analysis
[tissue_voxel_idx_list, nearest_cc_label] = fun_bin_data_to_idx_list(mask_labeled(tissue_nearest_mask_ind(valid_tissue_voxel_ind)));
clearvars tissue_nearest_mask_ind
% Quantify only the link specified by link_label: 
[process_link_label, list_idx_in_nearest_cc_label, list_idx_in_link_label] = intersect(nearest_cc_label, link_label);
tissue_voxel_idx_list = tissue_voxel_idx_list(list_idx_in_nearest_cc_label);
num_valid_label = numel(process_link_label);
% Count the number of voxels in each connected component
[mask_idx_list, mask_cc] = fun_bin_data_to_idx_list(mask_labeled);
mask_cc_size = cellfun(@numel, mask_idx_list) ./ (downsample_scale ^ im_dims);
[~, mask_cc_ind, tmp_link_label_idx] = intersect(mask_cc, link_label);
dtf.recon_mask_num_voxel(tmp_link_label_idx) = mask_cc_size(mask_cc_ind);
%% 
% tic_start = tic;
for iter_link_idx = 1 : num_valid_label
    tmp_link_list_idx = list_idx_in_link_label(iter_link_idx);
    tmp_link_label = link_label(tmp_link_list_idx);
    tmp_list_idx = tissue_voxel_idx_list{iter_link_idx};
    tmp_ind = valid_tissue_voxel_ind(tmp_list_idx);
    tissue_sub_g = fun_ind2sub(image_size, tmp_ind);
    % Has voxel on the boundary:
    tmp_on_boundaryQ = any(tissue_sub_g == 1, 'all') || any(tissue_sub_g == image_size, 'all');
    if ~tmp_on_boundaryQ
        tmp_dt = valid_tissue_space_dt_list(tmp_list_idx);
        dtf.tissue_dt_prctile(:, tmp_link_list_idx) = double(prctile(tmp_dt, tissue_dt_prctile));        
        dtf.nearest_tissue_volume(tmp_link_list_idx) = numel(tmp_dt) ./ (downsample_scale ^ im_dims);
        dtf.nearest_tissue_dt_mean(tmp_link_list_idx) = mean(tmp_dt);
        
        tmp_link_nearest_DT_max =  max(tmp_dt);
        dtf.nearest_tissue_dt_max(tmp_link_list_idx) = tmp_link_nearest_DT_max;
        dtf.nearest_tissue_dt_median(tmp_link_list_idx) = median(tmp_dt);
    %% Perturbation
        %% Compute distance transform perturbation
        % Determine the local bounding box that is sufficiently large for computing
        % the distance transform
        bbox_min = min(max(floor(tissue_sub_g - tmp_link_nearest_DT_max - bbox_expand_length), 1), [], 1);
        bbox_max = max(min(ceil(tissue_sub_g + tmp_link_nearest_DT_max + bbox_expand_length), image_size), [], 1);
        bbox_mmll = [bbox_min, bbox_max - bbox_min + 1];
        % Crop the label array and the distance transform array
        local_label_array = crop_bbox3(mask_labeled, bbox_mmll, 'default');
        % Remove the link from the local label array
        link_mask_ind_l = find(local_label_array == tmp_link_label);
        link_mask_sub_l = fun_ind2sub(bbox_mmll(4:6), link_mask_ind_l);
        if isempty(link_mask_ind_l)
            continue;
        end
        perturbed_mask = (local_label_array ~= 0);
        perturbed_mask(link_mask_ind_l) = false;
        [perturbed_dt, perturbed_nearest_mask_ind] = bwdist(perturbed_mask);
        % Rescale the distance
        perturbed_dt = perturbed_dt ./ downsample_scale;         
        
        tissue_sub_l = tissue_sub_g - bbox_min + 1;
        tissue_ind_l = sub2ind(bbox_mmll(4:6), tissue_sub_l(:, 1), ...
            tissue_sub_l(:, 2), tissue_sub_l(:, 3));
        % Combine the link indices and its nearest tissue indices
        % The DT properties of the rest of the voxel should not change, as they
        % are initially closer to other segment
        p_vol_ind_l = cat(1, tissue_ind_l, link_mask_ind_l);        
        tissue_vol_dt = perturbed_dt(tissue_ind_l);
        local_rm_link_dt = perturbed_dt(link_mask_ind_l);
        p_vol_dt = cat(1, tissue_vol_dt, local_rm_link_dt);
        %% Overall perturbation result        
        tissue_vol_dt_af_minus_bf = tissue_vol_dt - tmp_dt;
        dtf.diff_dp_d_mean(tmp_link_list_idx) = mean(tissue_vol_dt_af_minus_bf, 'omitnan');
        dtf.up_ts_2_p_ts_dt_mean(tmp_link_list_idx) = mean(tmp_dt ./ tissue_vol_dt, 'omitnan');
        %% New maximum in the perturbed space
        [p_vol_dt_max, p_vol_dt_max_list_ind] = max(p_vol_dt);
        % Position of the maximum DT in the perturbed space
        p_vol_dt_max_sub_l = fun_ind2sub(bbox_mmll(4:6), p_vol_ind_l(p_vol_dt_max_list_ind));
        % Distance between the removed link and the maximum DT voxel in the
        % perturbed space
        dist_rm_link_2_max_DT_in_pt_vol = pdist2(p_vol_dt_max_sub_l, link_mask_sub_l);
        min_dist_rm_link_2_max_DT = min(dist_rm_link_2_max_DT_in_pt_vol) ./ downsample_scale;
        
        dtf.pt_vol_dt_max(tmp_link_list_idx) = p_vol_dt_max;
        dtf.dist_rm_lk_2_pt_max_dt(tmp_link_list_idx) = min_dist_rm_link_2_max_DT;
        %% Average DT in the perturbed space
        dtf.pt_vol_dt_mean(tmp_link_list_idx) = mean(p_vol_dt);
        dtf.pt_vol_dt_median(tmp_link_list_idx) = median(p_vol_dt);        
        %% Maximum DT passed by the removed link after perturbation
        dtf.lk_dt_af_rm_mean(tmp_link_list_idx) = mean(local_rm_link_dt);
        dtf.lk_dt_af_rm_median(tmp_link_list_idx) = median(local_rm_link_dt);
        dtf.lk_dt_af_rm_max(tmp_link_list_idx) = max(local_rm_link_dt);
        %% Nearest Local maximum values and location with respect to the removed link
        is_local_max_Q_rm = fun_array_local_maximum(perturbed_dt, nonmax_win_size);
        p_vol_l_dt_max_ind = p_vol_ind_l(is_local_max_Q_rm(p_vol_ind_l));        
        if ~isempty(p_vol_l_dt_max_ind)
            roi_local_max_sub = fun_ind2sub(bbox_mmll(4:6), p_vol_l_dt_max_ind);
            dist_new_max_2_rm_link = pdist2(roi_local_max_sub, link_mask_sub_l);
            dist_new_max_2_rm_link_min = min(dist_new_max_2_rm_link, [], 2);
            [min_dist_2_local_max, nearest_local_max_idx] = min(dist_new_max_2_rm_link_min);
            % Distance between the removed link to the closest new local
            % maximum
            dtf.dist_rm_lk_2_nlm(tmp_link_list_idx) = min_dist_2_local_max ./ downsample_scale;
            % Value of the nearest local maximum after link removal
            dtf.nlm_v_af_rm(tmp_link_list_idx) = perturbed_dt(p_vol_l_dt_max_ind(nearest_local_max_idx));
        end
        %% Number of neighboring links and their assigned volume fraction
        % 1. Number of neighboring link segments
        % 2. Volume of the purturbed space assigned to each neighboring link
        % segment

        % Voxels in ROI are assigned to the remaining reconstructed mask voxel
        p_vol_new_nrst_mask_ind = perturbed_nearest_mask_ind(p_vol_ind_l);
        p_vol_new_nrst_link_label = local_label_array(p_vol_new_nrst_mask_ind);
        
        [p_nrst_link_vol_ind, nrst_link_label] = fun_bin_data_to_idx_list(p_vol_new_nrst_link_label);
        p_nrst_link_vol_size = cellfun(@numel, p_nrst_link_vol_ind);
        % Distribute the voxel belong to node to the rest of the connected links
        is_link_label_Q = (nrst_link_label > 0);
        if ~any(is_link_label_Q)
            continue;
        else
            for iter_label = 1 : numel(nrst_link_label)
                if ~is_link_label_Q(iter_label)
                    tmp_label = nrst_link_label(iter_label);
                    tmp_connected_link_label = vessel_graph.node.connected_link_label{- tmp_label};
                    tmp_connected_link_label = tmp_connected_link_label(tmp_connected_link_label ~= tmp_link_label);
                    [~, tmp1, ~] = intersect(nrst_link_label, tmp_connected_link_label);
                    p_nrst_link_vol_size(tmp1) = p_nrst_link_vol_size(tmp1) + p_nrst_link_vol_size(iter_label) ./ numel(tmp1);
                    p_nrst_link_vol_size(iter_label) = 0;
                end
            end
        end
        nrst_link_label = nrst_link_label(is_link_label_Q);
        p_nrst_link_vol_size = p_nrst_link_vol_size(is_link_label_Q);
        p_nrst_link_vol_ratio = p_nrst_link_vol_size ./ sum(p_nrst_link_vol_size);
        % Record
        dtf.num_nb_lk(tmp_link_list_idx) = numel(p_nrst_link_vol_ratio);
        dtf.nb_lk_vol_r_max(tmp_link_list_idx) = max(p_nrst_link_vol_ratio);
        dtf.nb_lk_vol_r{tmp_link_list_idx} = p_nrst_link_vol_ratio;
        dtf.nb_lk_label{tmp_link_list_idx} = nrst_link_label;
        %% Geodesic distance to the neighboring links
        % 3. In the connectivity graph, geodesic distance to the neighboring link
        %   a. Geodesic distance to the dominant neighboring links
        %   b. Geodesic distance to the farest link
        %   c. Geodesic distance to the link with largest re-assigned volume
        %   (ratio)
        num_nn_link = numel(nrst_link_label);        
        nn_link_node_lable = vessel_graph.link.connected_node_label(nrst_link_label, :);
        % Use one node for path-finding
        nn_link_node_label_path = sort(nn_link_node_lable, 2, 'descend');
        nn_link_node_label_unused = nn_link_node_label_path(:, 2);
        nn_link_node_label_path = nn_link_node_label_path(:, 1);
        rm_link_node_label = vessel_graph.link.connected_node_label(tmp_link_label, :);
        rm_link_node_label = rm_link_node_label(rm_link_node_label > 0 & rm_link_node_label <= max_valid_node_label);
        if ~any(rm_link_node_label)
            % Isolated link with two unconnected endpoints
            continue;
        end        
        % Get the label of the node of "branch order" 1
        ob1_node_label = cat(1, vessel_graph.node.connected_link_label{rm_link_node_label});
        ob1_node_label = vessel_graph.link.connected_node_label(ob1_node_label, :);
        ob1_node_label = setdiff(unique(ob1_node_label), rm_link_node_label);
        
        num_edges_to_nn_link = nan(num_nn_link, 1);
        tmp_num_source_node = numel(rm_link_node_label);
        tmp_source_node_label = rm_link_node_label(1);
        for iter_nn_link = 1 : num_nn_link
            tmp_target_node_label = nn_link_node_label_path(iter_nn_link);
            
            if any(nn_link_node_label_unused(iter_nn_link) == rm_link_node_label) && ...
                    any(tmp_target_node_label == ob1_node_label)                    
                num_edges_to_nn_link(iter_nn_link) = 0;
            elseif tmp_target_node_label == 0 || ...
                    tmp_target_node_label > max_valid_node_label % Probably node in the self-loop
                num_edges_to_nn_link(iter_nn_link) = nan;
            else
                [tmp_geodesic_path, tmp_geo_length] = shortestpath(connectivity_graph, ...
                    tmp_source_node_label, tmp_target_node_label, 'Method', 'unweighted');
                if tmp_num_source_node == 2
                    [tmp_geodesic_path_2, tmp_geo_length_2] = shortestpath(connectivity_graph, ...
                        rm_link_node_label(2), tmp_target_node_label, 'Method', 'unweighted');
                    if tmp_geo_length_2 < tmp_geo_length
                        tmp_geo_length = tmp_geo_length_2;
                        tmp_geodesic_path = tmp_geodesic_path_2;
                    end
                end                
                if any(tmp_geodesic_path == nn_link_node_label_unused(iter_nn_link))
                    % When the target node is the one (of two ends of a
                    % link) that is far away from the source node and
                    % therefore the path contained the other end of the
                    % node
                    tmp_geo_length = tmp_geo_length - 1;
                end
                num_edges_to_nn_link(iter_nn_link) = tmp_geo_length;
            end
        end
        nn_link_branching_order = num_edges_to_nn_link + 1;
        
        dtf.nb_lk_bch_od{tmp_link_list_idx} = nn_link_branching_order;
        
        [dtf.nb_lk_bch_od_max(tmp_link_list_idx), tmp_max_idx] = max(nn_link_branching_order);
        dtf.nb_lk_bch_od_max_vol_r(tmp_link_list_idx) = p_nrst_link_vol_ratio(tmp_max_idx);
        [~, tmp_max_idx] = max(p_nrst_link_vol_ratio);
        dtf.nb_lk_vol_r_max_bch_od(tmp_link_list_idx) = nn_link_branching_order(tmp_max_idx);
        dtf.nb_lk_bch_od_1_vol_r(tmp_link_list_idx) = sum(p_nrst_link_vol_ratio(nn_link_branching_order == 1));
    end
%     if mod(tmp_link_list_idx, 100) == 0
%         fprintf('Finish processing link %d/%d. Elapsed time is %f seconds\n', tmp_link_list_idx, num_valid_label, toc(tic_start));
%     end
end
dtf.tissue_dt_prctile = dtf.tissue_dt_prctile';
dtf = fun_struct2table(dtf);
end