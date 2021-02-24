set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'stack1';
grid_version = '512_cube_v1';
layer = 18;
block_idx_X = 7;
block_idx_Y = 7;
DataManager.visualize_block_location_in_itksnap(dataset_name, stack, 'cropped_enhanced', grid_version, layer, block_idx_X, block_idx_Y);
grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
vascBlock = DataManager.load_block_data(dataset_name, stack, grid_version, layer,block_idx_X,block_idx_Y);
avi_str = VideoWriter(DataManager.fp_constructor(dataset_name, stack, 'visualization', 'WholeBrain_stack1_512_cube_v1_18_7_7.avi',true));
avi_str.FrameRate = 10;
avi_str.Quality = 75;
open(avi_str);
for frame_idx = 1 : size(vascBlock, 3)
    writeVideo(avi_str, vascBlock(:,:,frame_idx));
end
close(avi_str)