DataManager = FileManager;
%% Determine the covering 1072 cubes
dataset_name = 'WholeBrain';
stack = 'ML20190124';
grid_version = '240_cube';
region_sub = [72, 512, 455];
region_sub_1um = region_sub .* 16;
grid_c = DataManager.load_grid(dataset_name, stack, '240_cube_combined_5_o_1');
grid_c_label = find(all(grid_c.bbox_xyz_mmll_pixel_list(:, 1:3) < region_sub_1um, 2) & ...
    all(grid_c.bbox_xyz_mmxx_pixel_list(:, 4:6) > region_sub_1um, 2));
grid_c_label = grid_c_label(1);
%%
grid_c_bbox_sub = grid_c.bbox_grid_sub_list(grid_c_label, :);
grid_c_bbox_mmxx_grid = grid_c.bbox_xyz_mmxx_grid_list(grid_c_label, :);

vessel_skel = DataManager.load_blocks_files('skel', dataset_name, stack, grid_version, ...
    grid_c_bbox_mmxx_grid(1):grid_c_bbox_mmxx_grid(4), ...
    grid_c_bbox_mmxx_grid(2):grid_c_bbox_mmxx_grid(5), ...
    grid_c_bbox_mmxx_grid(3):grid_c_bbox_mmxx_grid(6), 'single', '240_cube_rec');
%% Local the post-perfusion data 
dataset_name_iv = 'DKLab';
stack_iv = 'Rui20180622_perfusion';
vg_iv = DataManager.load_graph_in_block(dataset_name_iv, stack_iv, 'merged', 0, 0, 0);
%% Find the bounding box of the in vivo imaged volume in whole brain coordinate 
% Crop the skeleton and try rigid registration 
crop_center = [490 610 593];
crop_radius = 350;
crop_bbox_min = crop_center - [300, 350, 300];
crop_bbox_max = crop_center + [350, 350, 350];
crop_bbox_mmll = [crop_bbox_min, crop_bbox_max - crop_bbox_min + 1];
vessel_skel_cropped = crop_bbox3(vessel_skel, crop_bbox_mmll);

vg_pf = fun_skeleton_to_graph(vessel_skel_cropped);
vg_pf = fun_graph_add_radius(vg_pf, vessel_skel_cropped);

vg_pf_fp = './Demo/Data/Demo_vessel_graph.mat';
vg_iv_fp = './Demo/Data/Demo_post_perfusion_graph.mat';
DataManager.write_data(vg_pf_fp, vg_pf);
DataManager.write_data(vg_iv_fp, vg_iv);