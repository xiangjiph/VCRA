function output = fun_compute_block_data(grid_info, layer_list, option)


DataManager = FileManager;
if nargin < 2
    layer_list = grid_info.layer;
    option = struct;
end
if ~isfield(option, 'use_ssd')
    option.use_ssd = true;
end
    
num_layer = grid_info.num_grid_layer;
for layer_idx = layer_list
    ori_sec_idx1 = max( floor( (grid_info.bbox_z_mmxx{layer_idx}(1) - 1) * grid_info.downsample_rate) + 1, 1);
    ori_sec_idx2 = min( ori_sec_idx1 + (grid_info.bbox_z_mmll{layer_idx}(2) * grid_info.downsample_rate - 1) , grid_info.mask_info.ori_data_size(3));
    num_sec = ori_sec_idx2 - ori_sec_idx1 + 1;
    ori_sec_list = ori_sec_idx1 : ori_sec_idx2;
    file_to_use_list = cell(num_sec, 1);
    for tmp_file_id = 1 : num_sec
        file_to_use_list{tmp_file_id} = DataManager.fp_enhanced_image(dataset_name, stack,ori_sec_list(tmp_file_id) - 1, enhance_image_version);
    end
    if option.use_ssd
        fprintf('Copy %d files to SSD\n', numel(file_to_use_list));
        tic
        [tmpFolder, file_to_use_list] = DataManager.copy_to_scratch(file_to_use_list);
        toc
    end
    fprintf('Load data to RAM\n');% For non-downsampled sample
    tic
    image_stack = zeros([grid_1x.data_xy_size, grid_1x.bbox_z_mmll{layer_idx}(2)], 'uint8');
    for secInLayer = 1 : grid_1x.bbox_z_mmll{layer_idx}(2)
        image_stack(:,:,secInLayer) = DataManager.load_single_tiff(file_to_use_list{secInLayer});
    end
    if option.use_ssd
        DataManager.clear_scratch_tmp_folder(tmpFolder);
    end
    
    
    
end





end
