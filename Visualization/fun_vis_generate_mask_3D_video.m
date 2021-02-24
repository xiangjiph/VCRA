function fp = fun_vis_generate_mask_3D_video(mask, fp)



% vis_image = block_data;
% vis_mask = tmp_vessel_mask;
% fp = DataManager.fp_constructor(dataset_name, stack, 'visualization', 'test_visualization.avi',true);
disp('Generate 3D visualization');
tic
block_size = size(mask);
figure_size = 2 * max(block_size);
vis_mask_smooth = imclose(mask > 0, strel('sphere', 3));
fv = isosurface(vis_mask_smooth);
fv = reducepatch(fv, 0.1);
im = figure('Visible', 'off', 'Position', [0, 0, figure_size, figure_size], 'color', 'none');
ax = axes(im);
axis(ax, 'off', 'equal')
patch(ax, 'Vertices', fv.vertices, 'Faces', fv.faces, 'FaceColor', 'r', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.9, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
% xlabel('X/\mum', 'Interpreter', 'tex');
% ylabel('Y/\mum');
% zlabel('Z/\mum')
% axis(ax, 'off');
axis vis3d
camlight(ax)
clearvars fig_str
fprintf('Write vedio to %s\n', fp);
avi_str = VideoWriter(fp);
avi_str.FrameRate = 10;
avi_str.Quality = 90;
open(avi_str);
phi_list = 0 : -5 : -355;
num_phi = numel(phi_list);
for iter_phi = num_phi : -1 : 1 
    tmp_phi = phi_list(iter_phi);
%     [tmpX, tmpY, tmpZ] = pol2cart(tmp_phi*pi/180, 120, 120);
%     view(ax, [tmpX, tmpY, tmpZ]);
    view(ax, [tmp_phi + 30, 30]);
    tmp_fig = getframe(im);
    drawnow;
    writeVideo(avi_str, tmp_fig);
end
close(avi_str)
toc
end