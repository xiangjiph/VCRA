function cube_stat_str = fun_analysis_compute_local_reconstruction_statistics(wholebrain_stat_str)
% This implementation is too slow. Maybe because of broadcasting the entire
% cell feature array? 
% Number of cores to run the computation in parallel: 
num_core = 18;
%% Initialization
num_bbox = numel(wholebrain_stat_str.grid_ind);

[all_link_stat_cell, capillary_link_stat_cell, all_node_stat_cell] = deal(cell(num_bbox, 1));

[all_link_anisotropy_cell, all_link_ori_svd_1_cell, all_link_ori_vec_cell, ...
    capillary_link_isotropy_cell, capillary_ori_svd_cell, capillary_ori_vec_cell] = deal(cell(num_bbox, 1));
%%
% Create cell array for parfor. Distributing wholebrain_stat_str is too
% expansive. 
uni_ori_isotropy_info = load('./Metadata/uni_ori_isotropy.mat');
link_cell = wholebrain_stat_str.link_features;
node_cell = wholebrain_stat_str.node_features;
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
parfor iter_bbox = 1 : num_bbox
% for iter_bbox = 1 : num_bbox
    tmp_tic = tic;
    tmp_link_table = link_cell{iter_bbox};
    tmp_node_table = node_cell{iter_bbox};    
    if ~isempty(tmp_link_table)
        % Statistics for all the vessels
        all_link_stat_cell{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_link_table);
        % End-to-end vector - unweighted
        tmp_isotropy_str = fun_analysis_get_link_anisotropy(tmp_link_table.ep1_to_ep2_direction_vec, true);
        if ~isempty(tmp_isotropy_str.svd_min2max)
            tmp_uni_ori_mean = uni_ori_isotropy_info.min2max.mean(tmp_isotropy_str.num_data);
            tmp_uni_ori_std = uni_ori_isotropy_info.min2max.std(tmp_isotropy_str.num_data);
            
            tmp_uni_ori_svd_mean = uni_ori_isotropy_info.svd_1.mean(tmp_isotropy_str.num_data);
            tmp_uni_ori_svd_std = uni_ori_isotropy_info.svd_1.std(tmp_isotropy_str.num_data);
            
            all_link_anisotropy_cell{iter_bbox} = (tmp_isotropy_str.svd_min2max - tmp_uni_ori_mean) / tmp_uni_ori_std;
            all_link_ori_svd_1_cell{iter_bbox} = (tmp_isotropy_str.svd_value_ratio(1) - tmp_uni_ori_svd_mean) / tmp_uni_ori_svd_std;
            all_link_ori_vec_cell{iter_bbox} = tmp_isotropy_str.svd_max_vec;
        end
        
        % Statistics for capillary
        tmp_capillary_Q = ~tmp_link_table.is_large_vessel_Q;
        capillary_link_stat_cell{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_link_table, tmp_capillary_Q);
        
        tmp_isotropy_str = fun_analysis_get_link_anisotropy(tmp_link_table.ep1_to_ep2_direction_vec(tmp_capillary_Q, :), true);
        if ~isempty(tmp_isotropy_str.svd_min2max)
            tmp_uni_ori_mean = uni_ori_isotropy_info.min2max.mean(tmp_isotropy_str.num_data);
            tmp_uni_ori_std = uni_ori_isotropy_info.min2max.std(tmp_isotropy_str.num_data);
                        
            tmp_uni_ori_svd_mean = uni_ori_isotropy_info.svd_1.mean(tmp_isotropy_str.num_data);
            tmp_uni_ori_svd_std = uni_ori_isotropy_info.svd_1.std(tmp_isotropy_str.num_data);
            
            capillary_link_isotropy_cell{iter_bbox} = (tmp_isotropy_str.svd_min2max - tmp_uni_ori_mean) / tmp_uni_ori_std;
            capillary_ori_svd_cell{iter_bbox} = (tmp_isotropy_str.svd_value_ratio(1) - tmp_uni_ori_svd_mean) / tmp_uni_ori_svd_std;
            capillary_ori_vec_cell{iter_bbox} = tmp_isotropy_str.svd_max_vec;
        end
    end
    % Statistics for all the nodes
    if ~isempty(tmp_node_table)
        all_node_stat_cell{iter_bbox} = fun_analysis_get_graph_feature_basic_stat(tmp_node_table);
    end
%     fprintf('Finish processing cube %d. Elapsed time is %f seconds.\n', iter_bbox, toc(tmp_tic));
end
%%
cube_stat_str = struct;
cube_stat_str.all_link_stat = all_link_stat_cell;
cube_stat_str.capillary_stat = capillary_link_stat_cell;
cube_stat_str.all_node_stat = all_node_stat_cell;

cube_stat_str.all_link_isotropy = all_link_anisotropy_cell;
cube_stat_str.all_link_ori_svd_1 = all_link_ori_svd_1_cell;
cube_stat_str.all_link_ori_vec = all_link_ori_vec_cell;

cube_stat_str.capillary_isotropy = capillary_link_isotropy_cell;
cube_stat_str.capillary_ori_svd_1 = capillary_ori_svd_cell;
cube_stat_str.capillary_ori_vec = capillary_ori_vec_cell;
end