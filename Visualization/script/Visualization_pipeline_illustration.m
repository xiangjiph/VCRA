DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
grid_name = '240_cube';
recon_mask_ver = '240_cube_re';
grid_info = DataManager.load_grid(dataset_name, stack, grid_name);
output_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Paper');
%% Pick a 240-cube for visualizatin
% fun_vis_grid(dataset_name, stack, grid_name);
vis_cube_label = 28590;
vis_cube_grid_sub = grid_info.bbox_grid_sub_list(vis_cube_label, :);

vis_im = DataManager.load_block_data(dataset_name, stack, grid_name, vis_cube_grid_sub(1), ...
    vis_cube_grid_sub(2), vis_cube_grid_sub(3));
vis_mask = DataManager.load_block_mask(dataset_name, stack, grid_name, vis_cube_grid_sub(1), ...
    vis_cube_grid_sub(2), vis_cube_grid_sub(3));
vis_mask = fun_reconstruct_block_mask(vis_mask);
vis_mask_skel = bwskel(vis_mask);

vis_skel = DataManager.load_block_skl(dataset_name, stack, recon_mask_ver, vis_cube_grid_sub(1), ...
    vis_cube_grid_sub(2), vis_cube_grid_sub(3));
vis_mask_recon = fun_skeleton_reconstruction(vis_skel.ind, max(2, vis_skel.r), vis_skel.block_size);
vis_skel_array = false(vis_skel.block_size);
vis_skel_array(vis_skel.ind) = true;
vis_graph = fun_skeleton_to_graph(vis_skel.ind, vis_skel.block_size);

vis_itk_mask = uint8(vis_mask_recon);
vis_itk_mask(vis_graph.node.pos_ind) = 3;
vis_itk_mask(vis_graph.link.pos_ind) = 2;
DataManager.visualize_itksnap(vis_im, vis_itk_mask);
%% Image mask 
vis_mask_fv = isosurface(smooth3(vis_mask), false);
% vis_mask_fv = reducepatch(vis_mask_fv, 20000);
% fig_hdl = figure;
t = tiledlayout(1, 1);
% fig_hdl.Position(3:4) = 1.5 * fig_hdl.Position(3:4);
ax_hdl = nexttile;
patch_hdl = patch(vis_mask_fv, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 1, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
camlight(ax_hdl)
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
% ax_hdl.XLabel.String = 'X (\mum)';
% ax_hdl.YLabel.String = 'Y (\mum)';
% ax_hdl.ZLabel.String = 'Z (\mum)';
% ax_hdl.FontSize = 14;
ax_hdl.XAxis.Visible = 'off';
ax_hdl.YAxis.Visible = 'off';
ax_hdl.ZAxis.Visible = 'off';
hold(ax_hdl, 'on');
s_hdl = slice(ax_hdl, single(vis_im), 1, 1, 1);
s_hdl(1).EdgeColor = 'none';
s_hdl(2).EdgeColor = 'none';
s_hdl(3).EdgeColor = 'none';
ax_hdl.Colormap = gray;
t.Padding = 'none';
t.TileSpacing = 'none';
ax_hdl.Color = 'none';
fig_hdl = gcf;
tmp_fp = fullfile(output_folder, sprintf('Mask_240_cube_%d.png', vis_cube_label));
fun_print_image_in_several_formats(fig_hdl, tmp_fp, {'.png'});
% volumeViewer(vis_mask_recon)
%% Compare the skeleton before and after graph refinement and padding
im_colormap = gray(intmax(class(vis_im)));
vis_sec_list = 201 : 240;
vis_im_z_proj = max(vis_im(:, :, vis_sec_list), [], 3);
vis_mask_skel_z_proj = max(vis_mask_skel(:, :, vis_sec_list), [], 3);
vis_skel_array_z_proj = max(vis_skel_array(:, :, vis_sec_list), [], 3);
%% Image max projection
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
im_hdl = image(vis_im_z_proj, 'CDataMapping', 'direct');
ax_hdl.Colormap = im_colormap;
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.Visible = 'off';
tmp_fp = fullfile(output_folder, sprintf('Image_240_cube_%d_max_z_projection_sec_%d_to_%d.png', vis_cube_label, vis_sec_list(1), vis_sec_list(end)));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% Image
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.Visible = 'off';
imshowpair(vis_mask_skel_z_proj, vis_skel_array_z_proj);
tmp_fp = fullfile(output_folder, sprintf('Skel_ori_vs_final_240_cube_%d_max_z_projection_sec_%d_to_%d.png', vis_cube_label, vis_sec_list(1), vis_sec_list(end)));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% Overlay the skeleton im with the max-projection of the image
fig_hdl = figure;
ax_hdl1 = axes(fig_hdl);
im_hdl = imagesc(vis_im_z_proj, 'CDataMapping', 'direct');
ax_hdl1.DataAspectRatio = [1,1,1];
ax_hdl1.Visible = 'off';
ax_hdl1.Colormap = gray;
% Overlay with mask
ax_hdl3 = axes(fig_hdl);
tmp_vis_skel = zeros(size(vis_mask_skel_z_proj), 'uint16');
tmp_vis_skel(vis_mask_skel_z_proj & vis_skel_array_z_proj) = 1;
tmp_vis_skel(vis_mask_skel_z_proj & ~vis_skel_array_z_proj) = 2;
tmp_vis_skel(~vis_mask_skel_z_proj & vis_skel_array_z_proj) = 3;
im_hdl2 = image(ax_hdl3, tmp_vis_skel, 'CDataMapping', 'direct');
ax_hdl3.Visible = 'off';
linkprop([ax_hdl1, ax_hdl3], {'Position', 'XLim', 'YLim', 'DataAspectRatio', 'YDir', 'XAxisLocation', ...
    'YAxisLocation'});
im_hdl2.AlphaData = double(tmp_vis_skel > 0) .* 0.5;
ax_hdl3.Colormap = gray;
ax_hdl3.Colormap(2, :) = [1,0,0];
ax_hdl3.Colormap(3, :) = [0,0,1];
ax_hdl3.Colormap(4, :) = [0,1,0];
tmp_fp = fullfile(output_folder, sprintf('Im_w_skel_240_cube_%d_max_z_projection_sec_%d_to_%d.png', vis_cube_label, vis_sec_list(1), vis_sec_list(end)));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% Overlay mask and skeleton with image
tmp_fv = isosurface(vis_mask_recon, false);
skel_sub = fun_ind2sub(vis_skel.block_size, vis_skel.ind);

vis_im_312 = permute(vis_im, [3,1,2]);
vis_skel_array_312 = permute(vis_skel_array, [3,1,2]);
vis_mask_recon_312 = permute(vis_mask_recon, [3,1,2]);
vis_max_proj_sec = 1;
vis_max_sec = size(vis_im_312, 3);
if vis_max_proj_sec > 1
    tmp_fp = fullfile(output_folder, sprintf('Im_w_maks_n_skel_240_cube_%d_max_proj_%d.avi', vis_cube_label, vis_max_proj_sec));
else
    tmp_fp = fullfile(output_folder, sprintf('Im_w_maks_n_skel_240_cube_%d.avi', vis_cube_label));
end

avi_str = VideoWriter(tmp_fp);
avi_str.FrameRate = 10;
avi_str.Quality = 100;
open(avi_str);

for iter_inisec = 1 : vis_max_sec
    fprintf('Printing section %d\n', iter_inisec);
    vis_sec_list = iter_inisec : min(vis_max_proj_sec + iter_inisec - 1, vis_max_sec);
    vis_im_z_proj = max(vis_im_312(:, :, vis_sec_list), [], 3);
    vis_fin_recon = max(vis_mask_recon_312(:, :, vis_sec_list), [], 3);
    vis_skel_array_z_proj = max(vis_skel_array_312(:, :, vis_sec_list), [], 3);
    
    fig_hdl = figure('Visible', 'off');
    fig_hdl.Position(3:4) = [2048, 1024];
    ax_hdl0 = subplot(1,2,1);
    s_hdl = slice(ax_hdl0, single(vis_im), vis_sec_list(1), [], []);
    hold(ax_hdl0, 'on');
    patch_hdl = patch(tmp_fv, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
    scatter3(ax_hdl0, skel_sub(:, 2), skel_sub(:, 1), skel_sub(:, 3), 5, 'g', 'filled');
    ax_hdl0.Colormap = gray;
    s_hdl(1).EdgeColor = 'none';
    ax_hdl0.ZAxis.Direction = 'reverse';
    ax_hdl0.YAxis.Direction = 'reverse';
    ax_hdl0.DataAspectRatio = [1,1,1];
    ax_hdl0.XLim = [0,250];
    ax_hdl0.YLim = [0,250];
    ax_hdl0.ZLim = [0, 250];
    ax_hdl0.XLabel.String = 'X (\mum)';
    ax_hdl0.YLabel.String = 'Y (\mum)';
    ax_hdl0.ZLabel.String = 'Z (\mum)';
    ax_hdl0.FontSize = 12;
    ax_hdl0.FontWeight = 'bold';
    
%     ax_hdl1 = axes(fig_hdl);
    ax_hdl1 = subplot(1,2,2);
    im_hdl = imagesc(vis_im_z_proj, 'CDataMapping', 'direct');
    ax_hdl1.DataAspectRatio = [1,1,1];
    % ax_hdl1.Visible = 'off';
    ax_hdl1.XLabel.String = 'Y (\mum)';
    ax_hdl1.YLabel.String = 'Z (\mum)';
    ax_hdl1.FontSize = 12;
    ax_hdl1.FontWeight = 'bold';
    if vis_max_proj_sec > 1
        ax_hdl1.Title.String = sprintf('Max projection of section %d to %d', vis_sec_list(1), vis_sec_list(end));
    else
        ax_hdl1.Title.String = sprintf('Section %d', vis_sec_list);
    end
    ax_hdl1.Colormap = gray;
    % Overlay with mask
    ax_hdl3 = axes(fig_hdl);
    tmp_vis_skel = zeros(size(vis_mask_skel_z_proj), 'uint8');
    tmp_vis_skel(vis_fin_recon) = 1;
    tmp_vis_skel(vis_skel_array_z_proj) = 3;
    im_hdl2 = image(ax_hdl3, tmp_vis_skel, 'CDataMapping', 'direct');
    ax_hdl3.Visible = 'off';
    linkprop([ax_hdl1, ax_hdl3], {'Position', 'XLim', 'YLim', 'DataAspectRatio', 'YDir', 'XAxisLocation', ...
        'YAxisLocation'});
    im_hdl2.AlphaData = double(tmp_vis_skel > 0) .* 0.3;
    ax_hdl3.Colormap = gray;
    ax_hdl3.Colormap(2, :) = [1,0,0];
    ax_hdl3.Colormap(3, :) = [0,0,1];
    ax_hdl3.Colormap(4, :) = [0,1,0];    
    ax_hdl3.DataAspectRatio = [1,1,1];
    writeVideo(avi_str, getframe(fig_hdl));   
    delete(fig_hdl);
end
close(avi_str)
%% Visualize the reconstructed mask and the centerline 
tmp_fv = isosurface(smooth3(vis_mask_recon), false);
tmp_fv = reducepatch(tmp_fv, 10000);
fig_hdl = figure;
fig_hdl.Position(3:4) = 1.5 * fig_hdl.Position(3:4);
ax_hdl = axes(fig_hdl);
patch_hdl = patch(tmp_fv, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.35);
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
% ax_hdl.XLabel.String = 'X (\mum)';
% ax_hdl.YLabel.String = 'Y (\mum)';
% ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.XAxis.Visible = 'off';
ax_hdl.YAxis.Visible = 'off';
ax_hdl.ZAxis.Visible = 'off';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';

skel_sub = fun_ind2sub(vis_skel.block_size, vis_skel.ind);
hold(ax_hdl, 'on');
scatter3(ax_hdl, skel_sub(:, 2), skel_sub(:, 1), skel_sub(:, 3), 5, 'g', 'filled');
% ax_hdl.Color = [1,1,1];
s_hdl = slice(ax_hdl, single(vis_im), 1, 1, 1);
s_hdl(1).EdgeColor = 'none';
s_hdl(2).EdgeColor = 'none';
s_hdl(3).EdgeColor = 'none';
ax_hdl.Colormap = gray;
tmp_fp = fullfile(output_folder, sprintf('Recon_240_cube_%d.png', vis_cube_label));
fun_print_image_in_several_formats(fig_hdl, tmp_fp, {'.png'});
%% load brain mask and downsample it 
% whole_brain_mask = DataManager.load_data(fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_brain_d16x_annotated_mask.nii.gz'));
whole_brain_mask = DataManager.load_brain_mask(dataset_name, stack, 'whole_brain_d16x_annotated');
% whole_brain_im = DataManager.load_data(fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff'));
% DataManager.load_single_tiff('/data/Vessel/WholeBrain/ML_2018_08_15/processed_data/mask/whole_brain_d16x_annotated/WholeBrain_ML_2018_08_15_whole_brain_d16x_annotated_mask.tiff');
% volumeViewer(whole_brain_mask);
% Resize the brain 
wb_mask_ds_rate = 0.25;
wb_mask_rz = imresize3(uint8(whole_brain_mask), wb_mask_ds_rate, 'Method', 'nearest');
% Permute the axis 
wb_mask_rz = permute(wb_mask_rz, [3, 2, 1]);
wb_mask_rz = smooth3(wb_mask_rz);
wb_fv = isosurface(wb_mask_rz, false);
wb_fv = reducepatch(wb_fv);
%% Whole brain mask with 1 240-cube
x_tick_val_scaled = 0 : 2000 : 14000;
x_tick_val = x_tick_val_scaled ./ 16 .* wb_mask_ds_rate;

y_tick_val_scaled = 0 : 5000 : 10000;
y_tick_val = y_tick_val_scaled ./ 16 * wb_mask_ds_rate;

z_tick_val_scaled = 0 : 2000 : 12000;
z_tick_val = z_tick_val_scaled ./ 16 .* wb_mask_ds_rate;

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
patch_hdl = patch(wb_fv, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
ax_hdl.ZAxis.Direction = 'reverse';
ax_hdl.XAxis.Direction = 'normal';
ax_hdl.YLabel.String = 'Rostral - caudal (mm)';
ax_hdl.XLabel.String = 'Left - right (mm)';
ax_hdl.ZLabel.String = 'Ventral - dorsal (mm)';
ax_hdl.FontSize = 14;
view(ax_hdl, [60, 10]);

ax_hdl.XTick = x_tick_val;
ax_hdl.XTickLabel = arrayfun(@num2str, x_tick_val_scaled / 1000, 'UniformOutput', false);

ax_hdl.YTick = y_tick_val;
ax_hdl.YTickLabel = arrayfun(@num2str, y_tick_val_scaled / 1000, 'UniformOutput', false);

ax_hdl.ZTick = z_tick_val;
ax_hdl.ZTickLabel = arrayfun(@num2str, z_tick_val_scaled / 1000, 'UniformOutput', false);

tmp_fp = fullfile(output_folder, sprintf('Brain_mask.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);

test_block_mmxx = grid_info.bbox_xyz_mmxx_list(vis_cube_label, :); 
test_block_mmxx = test_block_mmxx([3 2 1 6 5 4]);
% test_block_mmxx = grid_info.bbox_xyz_mmxx_list; 
% Generate the polygons for patch 
test_block_mmxx = test_block_mmxx .* wb_mask_ds_rate ./ 16;
% test_block_mmxx = test_block_mmxx([2,3,1,5,6,4]);
[polyX, polyY, polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(test_block_mmxx);

hold(ax_hdl, 'on');
patch_hdl_2 = patch(ax_hdl, polyX, polyY, polyZ, 'r');
patch_hdl_2.FaceAlpha = 0.25;
  
tmp_fp = fullfile(output_folder, sprintf('Brain_mask_with_240_cube_%d.png', vis_cube_label));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% Whole brain mask with all 240-cube
test_block_mmxx = grid_info.bbox_xyz_mmxx_list; 
% Generate the polygons for patch 
test_block_mmxx = test_block_mmxx .* wb_mask_ds_rate ./ 16;
[polyX, polyY, polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(test_block_mmxx);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = axes(fig_hdl);
patch_hdl = patch(wb_fv, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
ax_hdl.ZAxis.Direction = 'reverse';
ax_hdl.XAxis.Direction = 'normal';
ax_hdl.Color = 'none';
ax_hdl.Visible = 'on';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';

ax_hdl.XTickLabel = num2str(str2double(ax_hdl.XTickLabel) .* 16 ./ wb_mask_ds_rate);
ax_hdl.YTickLabel = num2str(str2double(ax_hdl.YTickLabel) .* 16 ./ wb_mask_ds_rate);
ax_hdl.ZTickLabel = num2str(str2double(ax_hdl.ZTickLabel) .* 16 ./ wb_mask_ds_rate);

hold(ax_hdl, 'on');
patch_hdl_2 = patch(ax_hdl, polyX, polyY, polyZ, 'r');
patch_hdl_2.FaceAlpha = 0.25;
view(ax_hdl, [225, 30]);

tmp_fp = fullfile(output_folder, sprintf('Brain_mask_with_all_240_cube.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp, {'.png'});
%% Visualize the combined mask 2
combined_grid_version = '240_cube_combined_5_o_2';
grid_c = DataManager.load_grid(dataset_name, stack, combined_grid_version);
% fun_vis_grid(dataset_name, stack, combined_grid_version);
vis_combined_grid_label = 1601;
vis_combined_grid_sub = grid_c.bbox_grid_sub_list(vis_combined_grid_label, :);
vis_combined_grid_mmxx = grid_c.bbox_xyz_mmxx_pixel_list(vis_combined_grid_label, :);
vis_combined_grid_ori = vis_combined_grid_mmxx(1:3);
vis_combined_grid_ori = [vis_combined_grid_ori, vis_combined_grid_ori];

vis_combined_grid_mmxx = vis_combined_grid_mmxx - vis_combined_grid_ori;

vis_subgrid_mmxx_list = grid_c.internal_subgrid_bbox_mmxx{vis_combined_grid_sub(1), ...
    vis_combined_grid_sub(2), vis_combined_grid_sub(3)};
vis_subgrid_mmxx_list = vis_subgrid_mmxx_list - vis_combined_grid_ori;

[grid_polyX, grid_polyY, grid_polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(vis_combined_grid_mmxx);
[subgrid_polyX, subgrid_polyY, subgrid_polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(vis_subgrid_mmxx_list);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = axes(fig_hdl);
patch_hdl_1 = patch(ax_hdl, grid_polyX, grid_polyY, grid_polyZ, 'r');
patch_hdl_1.FaceAlpha = 0.25;
patch_hdl_1.LineWidth = 2;
hold(ax_hdl, 'on');
patch_hdl_2 = patch(ax_hdl, subgrid_polyX, subgrid_polyY, subgrid_polyZ, 'g');
patch_hdl_2.FaceAlpha = 0.25;
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.Color = 'none';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
ax_hdl.ZDir = 'normal';
view(ax_hdl, 3);
tmp_fp = fullfile(output_folder, sprintf('Grid_240_c5_o2_with_internal_grid.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% Visualize the combined mask 1 with all its internal grid
grid_c_1 = DataManager.load_grid(dataset_name, stack, '240_cube_combined_5_o_1');

vis_combined_grid_label = 800;
vis_combined_grid_sub = grid_c_1.bbox_grid_sub_list(vis_combined_grid_label, :);
vis_combined_grid_mmxx = grid_c_1.bbox_xyz_mmxx_pixel_list(vis_combined_grid_label, :);

vis_combined_grid_ori = vis_combined_grid_mmxx(1:3);
vis_combined_grid_ori = [vis_combined_grid_ori, vis_combined_grid_ori];

vis_combined_grid_mmxx = vis_combined_grid_mmxx - vis_combined_grid_ori;

vis_subgrid_mmxx_list = grid_c_1.sub_grid_bbox_mmxx{vis_combined_grid_sub(1), ...
    vis_combined_grid_sub(2), vis_combined_grid_sub(3)};
vis_subgrid_mmxx_list = vis_subgrid_mmxx_list - vis_combined_grid_ori;

[grid_polyX, grid_polyY, grid_polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(vis_combined_grid_mmxx);
[subgrid_polyX, subgrid_polyY, subgrid_polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(vis_subgrid_mmxx_list);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = axes(fig_hdl);
patch_hdl_1 = patch(ax_hdl, grid_polyX, grid_polyY, grid_polyZ, 'r');
patch_hdl_1.FaceAlpha = 0.0;
patch_hdl_1.LineWidth = 2;
hold(ax_hdl, 'on');
patch_hdl_2 = patch(ax_hdl, subgrid_polyX, subgrid_polyY, subgrid_polyZ, 'g');
patch_hdl_2.FaceAlpha = 0.4;
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.Color = 'none';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
ax_hdl.ZDir = 'normal';
view(ax_hdl, 3);

tmp_fp = fullfile(output_folder, sprintf('Grid_240_c5_o1_with_internal_grid.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% Brain mask with 1 c5_o1 grid
% Figure out which c5_o1 grid contains 240_cube vis_cube_label
vis_cube_bbox = grid_info.bbox_xyz_mmxx_list(vis_cube_label, :);
is_in_Q = all(grid_c_1.bbox_xyz_mmxx_pixel_list(:, 1:3) <= vis_cube_bbox(1:3), 2) & ...
    all(grid_c_1.bbox_xyz_mmxx_pixel_list(:, 4:6) >= vis_cube_bbox(4:6), 2);
vis_combined_grid_label = find(is_in_Q);

test_block_mmxx = grid_c_1.bbox_xyz_mmxx_pixel_list(vis_combined_grid_label, :);
% Generate the polygons for patch 
test_block_mmxx = test_block_mmxx .* wb_mask_ds_rate ./ 16;
[polyX, polyY, polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(test_block_mmxx);
[polyX_sub, polyY_sub, polyZ_sub] = fun_vis_bbox_mmxx_to_polygon_vertice(vis_cube_bbox .* wb_mask_ds_rate ./ 16);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = axes(fig_hdl);
patch_hdl = patch(wb_fv, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
ax_hdl.ZAxis.Direction = 'reverse';
ax_hdl.XAxis.Direction = 'normal';
ax_hdl.Color = 'none';
ax_hdl.Visible = 'on';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';

ax_hdl.XTickLabel = num2str(str2double(ax_hdl.XTickLabel) .* 16 ./ wb_mask_ds_rate);
ax_hdl.YTickLabel = num2str(str2double(ax_hdl.YTickLabel) .* 16 ./ wb_mask_ds_rate);
ax_hdl.ZTickLabel = num2str(str2double(ax_hdl.ZTickLabel) .* 16 ./ wb_mask_ds_rate);

hold(ax_hdl, 'on');
patch_hdl_2 = patch(ax_hdl, polyX, polyY, polyZ, 'g');
patch_hdl_2.FaceAlpha = 0.25;

patch_hdl_3 = patch(ax_hdl, polyX_sub, polyY_sub, polyZ_sub, 'r');
patch_hdl_3.FaceAlpha = 0.25;
view(ax_hdl, [225, 30]);

tmp_fp = fullfile(output_folder, sprintf('Brain_mask_with_240_cube_and_c5o1_grid.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp, {'.png'});
%% Visualize a layer of cubes
vis_cube_grid_layer = grid_info.bbox_grid_sub_list(vis_cube_label, 3);

vis_cube_bbox = grid_info.bbox_xyz_mmxx_list(vis_cube_label, :); 

test_block_mmxx = grid_info.bbox_xyz_mmxx{vis_cube_grid_layer};
% Generate the polygons for patch 
test_block_mmxx = test_block_mmxx .* wb_mask_ds_rate ./ 16;
[polyX, polyY, polyZ] = fun_vis_bbox_mmxx_to_polygon_vertice(test_block_mmxx);
[polyX_sub, polyY_sub, polyZ_sub] = fun_vis_bbox_mmxx_to_polygon_vertice(vis_cube_bbox .* wb_mask_ds_rate ./ 16);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = axes(fig_hdl);
patch_hdl = patch(wb_fv, 'FaceColor', [0.6, 0.6, 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.25);
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
ax_hdl.ZAxis.Direction = 'reverse';
ax_hdl.XAxis.Direction = 'normal';
ax_hdl.Color = 'none';
ax_hdl.Visible = 'on';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';

ax_hdl.XTickLabel = num2str(str2double(ax_hdl.XTickLabel) .* 16 ./ wb_mask_ds_rate);
ax_hdl.YTickLabel = num2str(str2double(ax_hdl.YTickLabel) .* 16 ./ wb_mask_ds_rate);
ax_hdl.ZTickLabel = num2str(str2double(ax_hdl.ZTickLabel) .* 16 ./ wb_mask_ds_rate);

hold(ax_hdl, 'on');
patch_hdl_2 = patch(ax_hdl, polyX, polyY, polyZ, 'g');
patch_hdl_2.FaceAlpha = 0.25;

% patch_hdl_3 = patch(ax_hdl, polyX_sub, polyY_sub, polyZ_sub, 'r');
% patch_hdl_3.FaceAlpha = 0.25;
view(ax_hdl, [225, 30]);

tmp_fp = fullfile(output_folder, sprintf('Brain_mask_with_a_layer_of_240_cubes.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp, {'.png'});
%% Annotation figures
fig_hdl = gcf;
tmp_fp = fullfile(output_folder, sprintf('Annotation_link_with_1_ep_1.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% 
fig_hdl = gcf;
% tmp_fp = fullfile(output_folder, sprintf('Annotation_link_dim_short_not_sure.png'));
% tmp_fp = fullfile(output_folder, sprintf('Annotation_link_dim_short_very_close5.png'));
tmp_fp = fullfile(output_folder, sprintf('Annotation_linker_wrong_1.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
%% After running Analysis_240_cube_whole_brain_stat_visualization_c5o2.m
%% Max projection of the reconstructed vessels in layer 37 on coronal plane
output_fp = fullfile(output_folder, sprintf('%s_%s_vessel_recon_maks_coronal_max_proj_layer_37.png', ...
    dataset_name, stack));
fun_print_image_in_several_formats(fig_handle, output_fp);
%%
output_fp = fullfile(output_folder, sprintf('%s_%s_vessel_recon_mask_with_vessel_length_density_layer_37.png', ...
    dataset_name, stack));
fun_print_image_in_several_formats(fig_handle, output_fp);
%% Illustration of capillary branch order
vessel_graph = fun_skeleton_to_graph(vis_skel.ind, vis_skel.block_size);
vessel_graph = fun_graph_add_radius(vessel_graph, ...
    sparse(double(vis_skel.ind), 1, double(vis_skel.r), prod(vis_skel.block_size), 1), 0.25);
vessel_graph.link.features = fun_analysis_get_link_features_by_label(vessel_graph);
is_capillary_Q = vessel_graph.link.features.dt_median <= 4;
vessel_graph.link.features.capillary_branch_order = fun_analysis_get_capillary_order_to_labeled_vessels(...
    vessel_graph, is_capillary_Q, true, false);
% Reconstruct vessel mask by capillary branching order
link_order = vessel_graph.link.features.capillary_branch_order + 1;
link_order(isnan(link_order)) = 1;
link_order = repelem(link_order, vessel_graph.link.num_voxel_per_cc, 1);
vis_mask_recon_order = fun_skeleton_reconstruction_label_aprox(...
    vessel_graph.link.pos_ind, full(vessel_graph.radius(vessel_graph.link.pos_ind)), ...
    link_order, vessel_graph.num.mask_size);
vis_mask_recon_order(vis_mask_recon_order > 3) = 0;
vis_cc = bwconncomp(vis_mask_recon_order > 0);
vis_cc_size = cellfun(@numel, vis_cc.PixelIdxList);
% Keep the largest one
vis_kept_Q = vis_cc_size == max(vis_cc_size);
vis_mask_recon_order(cat(1, vis_cc.PixelIdxList{~vis_kept_Q})) = 0;

vis_mask_fv_1 = isosurface(vis_mask_recon_order == 1, false);
vis_mask_fv_2 = isosurface(vis_mask_recon_order == 2, false);
vis_mask_fv_3 = isosurface(vis_mask_recon_order == 3, false);

fig_hdl = figure;
% fig_hdl.Position(3:4) = 3 * fig_hdl.Position(3:4);
ax_hdl = axes(fig_hdl);
patch_hdl = patch(vis_mask_fv_1, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
camlight(ax_hdl)
hold(ax_hdl, 'on');
patch_hdl_2 = patch(vis_mask_fv_2, 'FaceColor', 'blue', 'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
hold(ax_hdl, 'on');
patch_hdl_3 = patch(vis_mask_fv_3, 'FaceColor', 'green', 'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.XLim = [0, 240];
ax_hdl.YLim = [0, 240];
ax_hdl.ZLim = [0, 240];
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
legend(ax_hdl, [patch_hdl, patch_hdl_2, patch_hdl_3], 'Noncapillary', 'Branch order 1', 'Branch order 2', ...
    'Location', 'northeast')
tmp_fp = fullfile(output_folder, sprintf('Mask_240_cube_%d_example_capillary_branch_order.png', vis_cube_label));
fun_print_image_in_several_formats(fig_hdl, tmp_fp, {'.png'});
%% Illustrate distance transform
vis_dt_proj_sec = 80 : 160;
vis_dt_im = vis_im(:, :, vis_dt_proj_sec);
vis_dt_mask = vis_mask(:, :, vis_dt_proj_sec);
vis_dt_im_proj = max(vis_dt_im, [], 3);
vis_dt_mask_proj = max(vis_dt_mask, [], 3);
vis_dt_mask_dt = bwdist(vis_dt_mask_proj);
vis_dt_bbox_mmll = [67, 116, 60, 60];
vis_dt_crop_im = crop_bbox2(vis_dt_im_proj, vis_dt_bbox_mmll);
vis_dt_crop_dt = crop_bbox2(vis_dt_mask_dt, vis_dt_bbox_mmll);
vis_dt_crop_mask = crop_bbox2(vis_dt_mask_proj, vis_dt_bbox_mmll);
    %%
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
im_hdl = imagesc(ax_hdl, vis_dt_crop_im);
ax_hdl.DataAspectRatio = [1,1,1];
cbar_hdl = colorbar(ax_hdl);
cbar_hdl.Label.String = 'Intensity';
ax_hdl.FontSize = 14;
ax_hdl.XAxis.Visible = 'off';
ax_hdl.YAxis.Visible = 'off';
tmp_fp = fullfile(output_folder, sprintf('Illustration_dt_image.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);

fig_hdl = figure;
% fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1];
ax_hdl = axes(fig_hdl);
im_hdl = image(ax_hdl, vis_dt_crop_mask);
ax_hdl.Colormap = gray(2);
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.FontSize = 14;
ax_hdl.XAxis.Visible = 'off';
ax_hdl.YAxis.Visible = 'off';
tmp_fp = fullfile(output_folder, sprintf('Illustration_dt_mask.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);

fig_hdl = figure;
ax_hdl3 = axes(fig_hdl);
im_hdl_3 = imagesc(ax_hdl3, vis_dt_crop_dt);
ax_hdl3.DataAspectRatio = [1,1,1];
ax_hdl3.XTick = [0, 30, 60];
ax_hdl3.XAxis.Visible = 'off';
ax_hdl3.YAxis.Visible = 'off';
cbar_hdl_3 = colorbar(ax_hdl3);
cbar_hdl_3.Label.String = 'd (\mum)';
ax_hdl3.Colormap = jet;
ax_hdl3.FontSize = 14;
tmp_fp = fullfile(output_folder, sprintf('Illustration_dt_value.png'));
fun_print_image_in_several_formats(fig_hdl, tmp_fp);


% hold(ax_hdl2, 'on');
% im_hdl_3 = imagesc(ax_hdl2, vis_dt_crop_mask);
% im_hdl_3.CData = uint8(vis_dt_crop_mask) * 155;
% im_hdl_3.AlphaData = double(vis_dt_crop_mask) * 0.3;
%% Visualize a typical branch
% vis_node_label = 85;
% vis_node_ind = vis_graph.node.cc_ind{vis_node_label};
% vis_node_sub = fun_ind2sub(vis_graph.num.mask_size, vis_node_ind);
vis_node_sub = [90 61 132];
vis_node_ind = sub2ind(vis_graph.num.mask_size, vis_node_sub(1), vis_node_sub(2), vis_node_sub(3));
vis_node_label = full(vis_graph.node.map_ind_2_label(vis_node_ind));
vis_link_label = vis_graph.node.connected_link_label{vis_node_label};
vis_link_num_voxel = vis_graph.link.num_voxel_per_cc(vis_link_label);
vis_link_connected_node_label = vis_graph.link.connected_node_label(vis_link_label, :);
vis_link_connected_node_ind = cat(1, vis_graph.node.cc_ind{unique(vis_link_connected_node_label)});
vis_link_connected_node_sub = fun_ind2sub(vis_graph.num.mask_size, vis_link_connected_node_ind);
vis_bbox_mmxx = [min(vis_link_connected_node_sub, [], 1), max(vis_link_connected_node_sub, [], 1)];
vis_bbox_exp_half_length = 25;
vis_bbox_mmll = [vis_bbox_mmxx(1:3) - vis_bbox_exp_half_length, vis_bbox_mmxx(4:6) - vis_bbox_mmxx(1:3) + 1 + 2 * vis_bbox_exp_half_length];
% fun_vis_voxel_reconstruction_local(
vis_branch_mask = crop_bbox3(vis_mask_recon, vis_bbox_mmll);

vis_branch_mask_fv = isosurface(vis_branch_mask, false);
vis_branch_mask_fv = reducepatch(vis_branch_mask_fv);

t = tiledlayout(1, 1);
ax_hdl = nexttile;
patch_hdl = patch(vis_branch_mask_fv, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 1, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
camlight(ax_hdl)
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
ax_hdl.ZDir = 'reverse';
% ax_hdl.XLabel.String = 'X (\mum)';
% ax_hdl.YLabel.String = 'Y (\mum)';
% ax_hdl.ZLabel.String = 'Z (\mum)';
% ax_hdl.FontSize = 14;
ax_hdl.XAxis.Visible = 'off';
ax_hdl.YAxis.Visible = 'off';
ax_hdl.ZAxis.Visible = 'off';
t.Padding = 'none';
t.TileSpacing = 'none';
ax_hdl.Color = 'none';
ax_hdl.View = [90 -40];
fig_hdl = gcf;
fig_fp = fullfile(output_folder, 'Illustration_typical_branch_geometry.png');
fun_print_image_in_several_formats(fig_hdl, fig_fp);
