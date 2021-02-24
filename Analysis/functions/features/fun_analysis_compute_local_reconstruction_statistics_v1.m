function cube_stat_str = fun_analysis_compute_local_reconstruction_statistics(wholebrain_stat_str)

wb_grid_size = size(wholebrain_stat_str.volume_ratio);
%% Maximum radius of capillary
capillary_max_r = 3;
% Number of cores to run the computation in parallel: 
num_core = 18;
%% Initialization
[all_link_stat_cell, capillary_link_stat_cell, all_node_stat_cell]= deal(cell(wb_grid_size));
[all_link_isotropy_cell, all_link_ori_svd_1_cell, all_link_ori_vec_cell, ...
    capillary_link_isotropy_cell, capillary_ori_svd_cell, capillary_ori_vec_cell] = deal(cell(wb_grid_size(3), 1));
%%
% Create cell array for parfor. Distributing wholebrain_stat_str is too
% expansive. 
uni_ori_isotropy_info = load('./Metadata/uni_ori_isotropy.mat');
link_cell = wholebrain_stat_str.link_features;
node_cell = wholebrain_stat_str.node_features;
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
parfor iter_layer = 1 : wb_grid_size(3)
    tmp_tic = tic;
    fprintf('Processing layer %d/%d\n', iter_layer, wb_grid_size(3));
    % Initialization data for each layer for parfor
    layer_link_cell = link_cell(:, :, iter_layer);
    layer_node_cell = node_cell(:, :, iter_layer);
    [layer_all_link_stat_cell, layer_capillary_link_stat_cell, layer_node_stat_cell] = deal(cell(size(layer_link_cell)));
    [layer_all_link_anisotropy, layer_capillary_link_isotropy,...
        layer_all_link_ori_svd1, layer_capillary_ori_svd1] = deal(nan(size(layer_link_cell)));
    [layer_all_link_ori_vec, layer_capillary_ori_vec] = deal(nan([3, size(layer_link_cell)]));
    tmp_valid_ind = find(~cellfun(@isempty, layer_link_cell));
    num_valid_cell = numel(tmp_valid_ind);
    
    for iter_block = 1 : num_valid_cell
        tmp_ind = tmp_valid_ind(iter_block);
        
        tmp_link_table = layer_link_cell{tmp_ind};
        tmp_node_table = layer_node_cell{tmp_ind};
        if ~isempty(tmp_link_table)
            % Statistics for all the vessels
            layer_all_link_stat_cell{tmp_ind} = fun_analysis_get_graph_feature_basic_stat(tmp_link_table);
            tmp_isotropy_str = fun_analysis_get_link_anisotropy(tmp_link_table.ep1_to_ep2_direction_vec, true);
            if ~isempty(tmp_isotropy_str.svd_min2max)
                tmp_uni_ori_mean = uni_ori_isotropy_info.min2max.mean(tmp_isotropy_str.num_data);
                tmp_uni_ori_std = uni_ori_isotropy_info.min2max.std(tmp_isotropy_str.num_data);
                layer_all_link_anisotropy(tmp_ind) = (tmp_isotropy_str.svd_min2max - tmp_uni_ori_mean) / tmp_uni_ori_std;
                
                tmp_uni_ori_svd_mean = uni_ori_isotropy_info.svd_1.mean(tmp_isotropy_str.num_data);
                tmp_uni_ori_svd_std = uni_ori_isotropy_info.svd_1.std(tmp_isotropy_str.num_data);
                layer_all_link_ori_svd1(tmp_ind) = (tmp_isotropy_str.svd_value_ratio(1) - tmp_uni_ori_svd_mean) / tmp_uni_ori_svd_std;
                layer_all_link_ori_vec(:, tmp_ind) = tmp_isotropy_str.svd_max_vec;                
            end
            % Statistics for capillary
            tmp_capillary_Q = tmp_link_table.dt_median <= capillary_max_r;
            layer_capillary_link_stat_cell{tmp_ind} = fun_analysis_get_graph_feature_basic_stat(tmp_link_table, tmp_capillary_Q);
            tmp_isotropy_str = fun_analysis_get_link_anisotropy(tmp_link_table.ep1_to_ep2_direction_vec(tmp_capillary_Q, :), true);
            if ~isempty(tmp_isotropy_str.svd_min2max)
                tmp_uni_ori_mean = uni_ori_isotropy_info.min2max.mean(tmp_isotropy_str.num_data);
                tmp_uni_ori_std = uni_ori_isotropy_info.min2max.std(tmp_isotropy_str.num_data);
                layer_capillary_link_isotropy(tmp_ind) = (tmp_isotropy_str.svd_min2max - tmp_uni_ori_mean) / tmp_uni_ori_std;
                
                tmp_uni_ori_svd_mean = uni_ori_isotropy_info.svd_1.mean(tmp_isotropy_str.num_data);
                tmp_uni_ori_svd_std = uni_ori_isotropy_info.svd_1.std(tmp_isotropy_str.num_data);
                layer_capillary_ori_svd1(tmp_ind) = (tmp_isotropy_str.svd_value_ratio(1) - tmp_uni_ori_svd_mean) / tmp_uni_ori_svd_std;
                layer_capillary_ori_vec(:, tmp_ind) = tmp_isotropy_str.svd_max_vec;
            end
        end
        % Statistics for all the nodes
        if ~isempty(tmp_node_table)
            layer_node_stat_cell{tmp_ind} = fun_analysis_get_graph_feature_basic_stat(tmp_node_table);
        end
    end
    all_link_stat_cell(:, :, iter_layer) = layer_all_link_stat_cell;
    capillary_link_stat_cell(:, :, iter_layer) = layer_capillary_link_stat_cell;
    all_node_stat_cell(:, :, iter_layer) = layer_node_stat_cell;
    
    all_link_isotropy_cell{iter_layer} = layer_all_link_anisotropy;
    all_link_ori_svd_1_cell{iter_layer} = layer_all_link_ori_svd1;
    all_link_ori_vec_cell{iter_layer} = layer_all_link_ori_vec;
    
    capillary_link_isotropy_cell{iter_layer} = layer_capillary_link_isotropy;
    capillary_ori_svd_cell{iter_layer} = layer_capillary_ori_svd1;
    capillary_ori_vec_cell{iter_layer} = layer_capillary_ori_vec;
    
    fprintf('Finish processing layer %d. Elapsed time is %f\n', iter_layer, toc(tmp_tic));
end
%%
cube_stat_str = struct;
cube_stat_str.all_link_stat = all_link_stat_cell;
cube_stat_str.capillary_stat = capillary_link_stat_cell;
cube_stat_str.all_node_stat = all_node_stat_cell;

cube_stat_str.all_link_isotropy = cat(3, all_link_isotropy_cell{:});
cube_stat_str.all_link_ori_svd_1 = cat(3, all_link_ori_svd_1_cell{:});
cube_stat_str.all_link_ori_vec = cat(4, all_link_ori_vec_cell{:});

cube_stat_str.capillary_isotropy = cat(3, capillary_link_isotropy_cell{:});
cube_stat_str.capillary_ori_svd_1 = cat(3, capillary_ori_svd_cell{:});
cube_stat_str.capillary_ori_vec = cat(4, capillary_ori_vec_cell{:});
end