function wb_feature_skel_str = fun_analysis_load_whole_brain_240_cube_features_N_skel(grid_info, feature_version, skel_version, layer_list)

%% Initialization
persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
if nargin < 4
    layer_list = grid_info.layer;
end
num_core = 8;
same_file_Q = strcmp(feature_version, skel_version);
%% Simplify the parfor loading
load_bbox_grid_ind = cat(1, grid_info.bbox_grid_ind{layer_list});
load_bbox_grid_sub = cat(1, grid_info.bbox_grid_sub{layer_list});
num_load_bbox = size(load_bbox_grid_sub, 1);

[link_feature_cell, node_feature_cell,...
    skl_cell] = deal(cell(num_load_bbox, 1));
%%
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
parfor iter_bbox = 1 : num_load_bbox
    tmp_grid_sub = load_bbox_grid_sub(iter_bbox, :);
    % Load data
    try
        tmp_mask_str = DataManager.load_block_mask(dataset_name, stack, feature_version, ...
            tmp_grid_sub(1), tmp_grid_sub(2), tmp_grid_sub(3));
        try
            node_features = tmp_mask_str.node.features;
            link_features = tmp_mask_str.link.features;
            %% Further processing features here
            node_feature_cell{iter_bbox} = fun_analysis_postprocess_node_features(node_features, true);
            link_feature_cell{iter_bbox} = fun_analysis_postprocess_link_features(link_features);
            
            assert(isfield(tmp_mask_str, 'skl_ind'), 'The reconstruction file does not contain the skeleton data');
            % Need to check if more fields need to be removed to save memory
            if same_file_Q
                tmp_mask_str.ind = tmp_mask_str.skl_ind;
                tmp_mask_str = rmfield(tmp_mask_str, {'link', 'node', ...
                    'max_proj_1', 'max_proj_2', 'max_proj_3', 'stat', 'skl_ind'});
                skl_cell{iter_bbox} = tmp_mask_str;
            end
        catch ME1
            % For error handeling in parfor
            fprintf('Unable to process cube (%d, %d, %d)\n', tmp_grid_sub);
            rethrow(ME1);
        end
    catch ME
        fprintf('Unable to load feature in (%d, %d, %d)\n', tmp_grid_sub);
        fprintf('Error message: %s\n', ME.message);
        continue;
    end
    
    if ~same_file_Q
        try
            skl_cell{iter_bbox} = DataManager.load_block_skl(dataset_name, stack, skel_version, ...
                tmp_grid_sub(1), tmp_grid_sub(2), tmp_grid_sub(3));
        catch ME
            fprintf('Unable to load skeleton in (%d, %d, %d)\n', tmp_grid_sub);
            fprintf('Error message: %s\n', ME.message);
            continue;
        end
    end
end
pool_obj = gcp('nocreate');
delete(pool_obj);
%% Save to structure
wb_feature_skel_str = struct;
wb_feature_skel_str.grid_size = grid_info.grid_size;
wb_feature_skel_str.grid_ind = load_bbox_grid_ind;
wb_feature_skel_str.grid_sub = load_bbox_grid_sub;
wb_feature_skel_str.link_features = link_feature_cell;
wb_feature_skel_str.node_features = node_feature_cell;
wb_feature_skel_str.skl_str = skl_cell;
end