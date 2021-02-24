clc;clear; close all;
%% Graph analysis settings
opt_graph_feature = struct;
computeQ = struct;
computeQ.basic_link_feature = true;
computeQ.basic_node_feature = true;
computeQ.basic_reconstruction = true;
computeQ.node_path_to_nearest_neighbor = true;
computeQ.link_shortest_path = true;
computeQ.dist_tissue_2_vessel = true;
computeQ.dimension = false;
computeQ.max_z_proj = true;
computeQ.capillary_branching_order = false;
computeQ.capillary_artery_order = NaN; % Placeholder for future development
computeQ.capillary_vein_order = NaN; % Placeholder for future development

computeQ.noncapillary_DT = true; % Pick the vessels are are not capillaries, analyze the distance transform properties.
computeQ.sgl_seg_rm_ptb = true;
computeQ.dist_to_brain_surface = false; % Compute the distance between the vessels and the surface of the brain.

computeQ.ep2ep_anisotropy = false;
computeQ.volume_weighted_ep2ep_anisotropy = false;
computeQ.capillary_ep2ep_anisotropy = false;
computeQ.capillary_volume_weighted_ep2ep_anisotropy = false;

opt_graph_feature.computeQ = computeQ;
opt_graph_feature.nonmax_win_size_um = 20;
opt_graph_feature.vis_dim_fit = false;
opt_graph_feature.dim_fit_cutoff_length = 25;
opt_graph_feature.recon_max_error_rate = 0.1;
opt_graph_feature.merge_neighbor_nodes_Q = false;
opt_graph_feature.max_merge_length = 5;
opt_graph_feature.always_merge_max_length = 1;
opt_graph_feature.DT_scale_factor = 1;

opt_graph_feature.capillary_max_radius = 3.5; % For saperating capillary vs large vessels for local statistics.
% Overwrite the existing features in the graph (if exist)
opt_graph_feature.overwrite_computed_featureQ = false;
%% Grid indices
recon_r = 0.25;
lattice_name_list = {'(10, 3)-a', '(10, 3)-b', '(10, 3)-c', '(8, 3)-a', 'cubic'};
lattice_shortest_loop_num_edge_list = [10, 10, 10, 8, 4];
simu_bond_length_list = [30, 50, 70, 90];
simu_system_size_list = [8, 10];
num_lattice = numel(lattice_name_list);
num_bond_lenght = numel(simu_bond_length_list);
num_system_n = numel(simu_system_size_list);
[lattice_stat_cell, lattice_data_cell] = deal(cell(num_bond_lenght, num_system_n, num_lattice));
%%
task_tic = tic;
for iter_lattice = 1 : num_lattice
    % lattice_name = '(10, 3)-c';
    lattice_name = lattice_name_list{iter_lattice};
    lattice_parameters = fun_simulaton_get_lattice_parameters(lattice_name);
    for iter_system_n = 1 : num_system_n
        for iter_length = 1 : num_bond_lenght
            tmp_tic = tic;
            edge_length = simu_bond_length_list(iter_length);
            system_size = simu_system_size_list(iter_system_n);
            system_size = system_size * edge_length;
            %% Initialization
            tmp_stat_str = struct;
            tmp_stat_str.lattice_name = lattice_name;
            tmp_stat_str.set_length_pxl = edge_length;
            tmp_stat_str.system_size_0 = system_size;
            %% Parameters
            scale_factor = single(edge_length / lattice_parameters.bond_length);
            system_size_exp = system_size + sqrt(3) * scale_factor;
            n_in_max = (lattice_parameters.translation_vectors .* scale_factor)  \ [system_size; system_size; system_size];
            n_in_max = fix(n_in_max);
            n_expand_max = zeros(size(n_in_max));
            n_expand_min = zeros(size(n_in_max));
            n_vecs = cell(3, 1);
            for iter_dim = 1 : 3
                if n_in_max(iter_dim) > 0
                    n_expand_max(iter_dim) = n_in_max(iter_dim) + 3;
                    n_expand_min(iter_dim) = -3;
                else
                    n_expand_max(iter_dim) = 3;
                    n_expand_min(iter_dim) = n_in_max(iter_dim) - 3;
                end
                tmp_vec = n_expand_min(iter_dim) : n_expand_max(iter_dim);
                n_vecs{iter_dim} = tmp_vec;
            end
            % Compute site coordinate by translation
            [N1, N2, N3] = ndgrid(n_vecs{:});
            N_mat = cat(2, N1(:), N2(:), N3(:))';
            num_base_vec = size(lattice_parameters.base_vectors, 2);
            site_sub_cell = cell(num_base_vec, 1);
            for iter_vec = 1 : num_base_vec
                site_sub_cell{iter_vec} = (lattice_parameters.base_vectors(:, iter_vec) + lattice_parameters.translation_vectors * N_mat)';
            end
            site_sub = cat(1, site_sub_cell{:});
            %% Site selection
            site_sub_pxl = round(site_sub .* scale_factor + 1);
            tmp_in_Q = all(site_sub_pxl >= 1, 2) & all(site_sub_pxl <= system_size, 2);
            site_in_sub = site_sub(tmp_in_Q, :);
            site_in_sub_pxl = site_sub_pxl(tmp_in_Q, :);
            %% Determine the connections
            assert(size(site_sub, 2) == size(unique(site_sub, 'rows', 'stable'), 2));
            num_in_sites = size(site_in_sub, 1);
            site_pdist = pdist2(site_in_sub, site_in_sub);
            site_pdist = abs(site_pdist - lattice_parameters.bond_length) < 1e-2;
            [ind1, ind2] = find(site_pdist);
            connected_list_ind_pair = sort(cat(2, ind1, ind2), 2, 'ascend');
            connected_list_ind_pair = unique(connected_list_ind_pair, 'rows', 'stable');            
            %% Construct the lattice mask 
            lattice_mask = false(ones(1, 3) .* system_size);
            mask_size = size(lattice_mask);
            site_ind = sub2ind(mask_size, site_in_sub_pxl(:, 1), site_in_sub_pxl(:, 2), site_in_sub_pxl(:, 3));
            lattice_mask(site_ind) = true;
            connected_line_ind_cell = cell(size(connected_list_ind_pair));
            for iter_line = 1 : size(connected_line_ind_cell, 1)
                tmp_sub = site_in_sub_pxl(connected_list_ind_pair(iter_line, :), :);
                p2p_vec_ind = fun_graph_get_p2p_line_mask(mask_size, tmp_sub(1, :), tmp_sub(2, :));
                %         p2p_vec_ind = fun_graph_get_p2p_cylinder_mask(mask_size, tmp_sub(1, :), tmp_sub(2, :), 1);
                lattice_mask(p2p_vec_ind) = true;
            end
            tmp_stat_str.system_size = size(lattice_mask);
            %% Dilate the mask, skeletonization, recentering
            lattice_mask_d = imdilate(lattice_mask, strel('cube', 3));
            lattice_mask_d_dt = bwdist(~lattice_mask_d);
            lattice_mask_skl_d = bwskel(lattice_mask_d);
            lattice_mask_skl_d = fun_skeleton_recentering(lattice_mask_skl_d, lattice_mask_d_dt);
%             
%             vis_mask = uint8(lattice_mask_d);
%             vis_mask(lattice_graph.link.pos_ind) = 2;
%             vis_mask(lattice_graph.node.pos_ind) = 3;
%             DataManager.visualize_itksnap(uint8(lattice_mask_d), vis_mask);
            %% Skeletonization
%             lattice_mask_skel = bwskel(lattice_mask);
%             lattice_graph = fun_skeleton_to_graph(lattice_mask_skel);
            lattice_graph = fun_skeleton_to_graph(lattice_mask_skl_d);
            lattice_graph_ind = cat(1, lattice_graph.link.pos_ind, lattice_graph.node.pos_ind, ...
                lattice_graph.endpoint.pos_ind);
            lattice_graph.radius = sparse(lattice_graph_ind, 1, recon_r, lattice_graph.num.block_voxel, 1);
            
            [lattice_graph, recon_mask] = fun_analysis_compute_graph_features_by_labels(lattice_graph, [], [], opt_graph_feature);           
            %% Process analysis result            
            tmp_selected_link_Q = lattice_graph.link.features.shortest_loop_geodesic_length == lattice_shortest_loop_num_edge_list(iter_lattice) & ...
                isfinite(lattice_graph.link.features.nearest_tissue_dt_max);
            tmp_internal_bbox_mmxx = [ones(1, 3) *  edge_length / 2, ones(1, 3) * (system_size - edge_length / 2)];
            tmp_is_internal_link_Q = fun_voxel_sub_in_bbox_mmxx_Q(lattice_graph.link.features.ep1_sub, tmp_internal_bbox_mmxx) & ...
                fun_voxel_sub_in_bbox_mmxx_Q(lattice_graph.link.features.ep2_sub, tmp_internal_bbox_mmxx);
            
            tmp_record_feature_table = lattice_graph.link.features(tmp_selected_link_Q & tmp_is_internal_link_Q, :);
            tmp_is_outlier_Q = fun_analysis_is_outlier_by_percentile(tmp_record_feature_table.nearest_tissue_dt_max);
            tmp_record_feature_table = tmp_record_feature_table(~tmp_is_outlier_Q, :);
            tmp_record_feature_table = fun_getfields(tmp_record_feature_table, {'length', ...
                'ep2ep_dist', 'straightness', 'shortest_loop_length', 'shortest_loop_geodesic_length', ...
                'nearest_tissue_dt_max', 'nearest_tissue_dt_median', 'nearest_tissue_volume', 'nearest_tissue_dt_mean', ...
                'recon_mask_num_voxel', 'tissue_dt_prctile', ...
                'num_nb_lk', 'nb_lk_vol_r_max', 'nb_lk_bch_od_1_vol_r', ...
                'nb_lk_vol_r_max_bch_od', 'lk_dt_af_rm_mean', ...
                'lk_dt_af_rm_median', 'lk_dt_af_rm_max', 'pt_vol_dt_max',...
                'dist_rm_lk_2_pt_max_dt', 'pt_vol_dt_mean', 'pt_vol_dt_median', ...
                'up_ts_2_p_ts_dt_mean', 'diff_dp_d_mean' 'nearest_tissue_radius', ...
                'dist_rm_lk_2_nlm', 'nlm_v_af_rm'});            
            tmp_record_feature_table = struct2table(tmp_record_feature_table);
            tmp_record_feature_table.delta_dt_max = tmp_record_feature_table.pt_vol_dt_max - ...
                tmp_record_feature_table.nearest_tissue_dt_max;
            
            tmp_stat_str.recon_r = recon_r;
            tmp_stat_str.all_link_stat = fun_simulation_SFL_get_statistics(tmp_record_feature_table, tmp_stat_str.set_length_pxl);
                        
            tmp_is_nan_Q = any(isnan(table2array(tmp_record_feature_table)), 2);
            tmp_stat_str.dt_lm_link_stat = fun_simulation_SFL_get_statistics(tmp_record_feature_table(~tmp_is_nan_Q, :), tmp_stat_str.set_length_pxl);
            
%             tmp_stat_str.nearest_tissue_dt_max_mean * sqrt(tmp_stat_str.TargetTotalLengthDensity)
%             tmp_stat_str.nearest_tissue_dt_mean_mean * sqrt(tmp_stat_str.TargetTotalLengthDensity)
            %% Output
            fprintf('#######################################################\n');
            fprintf('Set bond length: %d\nSystem length: %d\n', edge_length, tmp_stat_str.system_size(1));
            fprintf('Finish this iteration. Elapsed time is %f seconds\n', toc(tmp_tic));
            fprintf('#######################################################\n');            
            %%
            lattice_stat_cell{iter_length, iter_system_n, iter_lattice} = tmp_stat_str;
            lattice_data_cell{iter_length, iter_system_n, iter_lattice} = lattice_graph;
        end
    end
end
%%
% Save result
DataManager = FileManager;
result_fp = DataManager.fp_analysis_data_file('WholeBrain', 'ML_2018_08_15', ...
    sprintf('Lattice_space_filling_link_properties_table_min_r_%.2f_um.mat', recon_r));
DataManager.write_data(result_fp, lattice_data_cell);
stat_fp = DataManager.fp_analysis_data_file('WholeBrain', 'ML_2018_08_15', ...
    sprintf('Lattice_space_filling_stat_table_min_r_%.2f_um.mat', recon_r));
DataManager.write_data(stat_fp, lattice_stat_cell);
fprintf('Finish task. Elapsed time is %f seconds.\n', toc(task_tic));