function proj_mask_str = fun_analysis_load_whole_brain_mask_max_proj(grid_info, recon_mask_name, layer_list)
% Currently only support loading data on the coronal (z) plane. 
% To do list: 
%   1. Allow loading max projectiong in any of the 3 directions from any
%   plane. 

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end
dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
if nargin < 3
    layer_list = grid_info.layer;
end
num_layer_to_process = numel(layer_list);
grid_sub_to_load = cat(1, grid_info.bbox_grid_sub{layer_list});
grid_sub_to_load = grid_sub_to_load.';
grid_ind_to_load = cat(1, grid_info.bbox_grid_ind{layer_list});
[grid_max_proj_1, grid_max_proj_2, grid_max_proj_3] = deal(cell(size(grid_ind_to_load)));
num_grid_cube = numel(grid_ind_to_load);
%% Parallel setting
num_core = 8;
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
%%
tmp_all_tic = tic;
parfor iter_cube = 1 : num_grid_cube
        try
            tmp_grid_sub = grid_sub_to_load(:, iter_cube);
            tmp_mask_str = DataManager.load_block_mask(dataset_name, stack, recon_mask_name, ...
                tmp_grid_sub(1), tmp_grid_sub(2), tmp_grid_sub(3));
            grid_max_proj_1{iter_cube} = tmp_mask_str.max_proj_1;
            grid_max_proj_2{iter_cube} = tmp_mask_str.max_proj_2;
            grid_max_proj_3{iter_cube} = tmp_mask_str.max_proj_3;
        catch ME
            fprintf('Unable to load reconstructed mask in grid (%d, %d, %d)\n', tmp_grid_sub);
            fprintf('Error message: %s\n', ME.message);
            continue;
        end
end
fprintf('Finish loading max projection of %s mask. Elapsed time is %f seconds\n', recon_mask_name, toc(tmp_all_tic));
%%
proj_mask_str = struct;
wb_grid_size = [grid_info.grid_size(1:2), num_layer_to_process];
[proj_mask_str.mask_max_proj_1, proj_mask_str.mask_max_proj_2, ...
    proj_mask_str.mask_max_proj_3] = deal(cell(wb_grid_size));

proj_mask_str.mask_max_proj_1(grid_ind_to_load) = grid_max_proj_1;
proj_mask_str.mask_max_proj_2(grid_ind_to_load) = grid_max_proj_2;
proj_mask_str.mask_max_proj_3(grid_ind_to_load) = grid_max_proj_3;

end