function [exitcode, varargout] = fun_segmentation_image_to_mask(dataset_name, stack, grid_version, image_grid_label, opt)

exitcode = -1; %#ok<NASGU>
DataManager = FileManager;
if isfield(opt, 'grid_info')
    grid_info = opt.grid_info;
else
    grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
end
% if ~isfield(opt, 'fun_handle')
%     seg_fun = @fun_mouselight_segmentation_1um_cube;
% else
%     seg_fun = opt.fun_handle;
% end
grid_idx_1 = grid_info.bbox_grid_sub_list(image_grid_label, 1);
grid_idx_2 = grid_info.bbox_grid_sub_list(image_grid_label, 2);
grid_layer = grid_info.bbox_grid_sub_list(image_grid_label, 3);
%% 
mask_str = struct;
mask_str.dataset_name = dataset_name;
mask_str.stack = stack;
mask_str.version = grid_version;
mask_str.global_bbox_mmll = grid_info.bbox_xyz_mmll_list(image_grid_label, :);
mask_str.global_bbxx_mmxx = grid_info.bbox_xyz_mmxx_list(image_grid_label, :);
mask_str.global_block_size = grid_info.data_size;
mask_str.grid_idx_1 = grid_idx_1;
mask_str.grid_idx_2 = grid_idx_2;
mask_str.grid_layer = grid_layer;
%% Processing
% Load data
block_data = DataManager.load_block_data(dataset_name, stack, grid_version, grid_idx_1, grid_idx_2, grid_layer);
% Generate segmentation
[vessel_mask, mask_str.record] = fun_mouselight_segmentation_1um_cube(block_data, opt);
% Convert mask to ind list
assert(numel(vessel_mask) < intmax('uint32'));
mask_str.ind = uint32(find(vessel_mask));
mask_str.block_size = size(vessel_mask);
% Save result
DataManager.write_block_mask(mask_str, dataset_name, stack, grid_version, ...
    grid_idx_1, grid_idx_2, grid_layer);
exitcode = 1;
if nargout == 2
    varargout{1} = mask_str;
end
end