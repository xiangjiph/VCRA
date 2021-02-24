function varargout = fun_vis_single_section_local_network_orientation_by_color_code(im, ori_vec_3d_list, patch_bbox_mmxx, vis_opt)
% To be implemented
% 1. Check the position of the colorbar. Overlap between the image and the
% overlap


%%
% Check if the patch value is nan:
is_nan_Q = isnan(ori_vec_3d_list);
is_nan_Q = all(is_nan_Q, 2);
num_valid_patch = nnz(~is_nan_Q);
% Remvoe nan patch
ori_vec_3d_list = ori_vec_3d_list(~is_nan_Q, :);
% Compute bounding box patches vertices
bbox_vertices_X = patch_bbox_mmxx(~is_nan_Q, [2,4,4,2])';
bbox_vertices_Y = patch_bbox_mmxx(~is_nan_Q, [1, 1, 3, 3])';
%% Pixel size
if isfield(vis_opt, 'pixel_size') && isscalar(vis_opt.pixel_size)
    pixel_size = vis_opt.pixel_size;
else
    pixel_size = 1;
end
%% Colorbar and color map
patch_v = abs(ori_vec_3d_list);
if isfield(vis_opt, 'vis_figQ') && islogical(vis_opt.vis_figQ)
    vis_fig_Q = vis_opt.vis_figQ;
else
    vis_fig_Q = false;
end
% im_size = size(im);
im_size = [2700, 3600];
if vis_fig_Q
    fig_handle = figure('Position', [1, 1, im_size(2), im_size(1)], 'Visible', 'on');
else
    fig_handle = figure('Position', [1, 1, im_size(2), im_size(1)], 'Visible', 'off');
end
% Set up the axes object for overlay
im_ax = axes(fig_handle);
image(im_ax, im, 'CDataMapping', 'scaled');
im_ax.Colormap = gray;
im_ax.Position = [0.075, 0.1, 0.8, 0.8];
im_ax.XLabel.String = 'Length/\mum';
im_ax.YLabel.String = 'Length/\mum';
im_ax.XAxis.TickLabels = cellfun(@(x) num2str(x, '%d'), ...
    num2cell(im_ax.XAxis.TickValues .* pixel_size), 'UniformOutput', false);
im_ax.YAxis.TickLabels = cellfun(@(x) num2str(x, '%d'), ...
    num2cell(im_ax.YAxis.TickValues .* pixel_size), 'UniformOutput', false);
im_ax.DataAspectRatio = [1, 1, 1];
ax1 = axes(fig_handle);
linkprop([im_ax, ax1], {'Position', 'XLim', 'YLim', 'DataAspectRatio', 'YDir', 'XAxisLocation', ...
    'YAxisLocation'});
ax1.Visible = 'off';
%% Draw patches
if num_valid_patch > 0
    p = patch(ax1, bbox_vertices_X, bbox_vertices_Y, ones(num_valid_patch,1), 'FaceAlpha', 0.5);
    % Directly map the direction to the 
    p.FaceVertexCData = patch_v;
    p.FaceAlpha = 0.6;
    p.CDataMapping = 'direct';
    p.FaceColor = 'flat';
end
%% Axes setting
if isfield(vis_opt, 'fig_title') && ischar(vis_opt.fig_title)
    im_ax.Title.String = vis_opt.fig_title;
    im_ax.FontSize = 24;
end
set(fig_handle, 'PaperPositionMode', 'auto');
%% Write image to folder
if isfield(vis_opt, 'file_path') && ~isempty(vis_opt.file_path)
    if isfield(vis_opt, 'resolution') && isscalar(vis_opt.resolution)
        fig_resolution = sprintf('-r%d', vis_opt.resolution);
    else
        fig_resolution = '-r300';
    end   
    print(fig_handle, '-dpng', '-opengl', fig_resolution, vis_opt.file_path);
end
if nargout == 1 
    if isfield(vis_opt, 'return_frameQ') && vis_opt.return_frameQ
        varargout{1} = getframe(fig_handle);
        close(fig_handle);
    elseif vis_fig_Q
        varargout{1} = fig_handle;
    end
elseif ~vis_fig_Q
    close(fig_handle);
    varargout{1} = 0;
end

end