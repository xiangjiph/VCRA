function [ax_hdl, varargout] = fun_vis_labeled_array(labeled_array, colormap_array, fig_hdl)


if nargin < 3
    fig_hdl = figure;
end
% Gather the indices of the labeled array
unique_label = unique(labeled_array);
is_vis_label_Q = unique_label ~= 0;
unique_label = unique_label(is_vis_label_Q);
num_vis_label = numel(unique_label);
fv_cell = cell(num_vis_label, 1);
tmp_tic = tic;
for iter_label = 1 : num_vis_label
    tmp_label_mask = smooth3(labeled_array == unique_label(iter_label));
    fv_cell{iter_label} = isosurface(tmp_label_mask, false);
    fv_cell{iter_label} = reducepatch(fv_cell{iter_label}, 0.1);
end
fprintf('Finish generating isosurface patches. Elapse time was %f seconds.\n', ...
    toc(tmp_tic));

if nargin < 2
    patch_coloer_map = lines(num_vis_label);
else
    patch_coloer_map = colormap_array;
end

ax_hdl = axes(fig_hdl);
for iter_label = 1 : num_vis_label
    tmp_hdl = patch(ax_hdl, fv_cell{iter_label}, 'FaceColor', patch_coloer_map(iter_label, :), ...
        'EdgeColor', 'none', 'FaceAlpha', 0.9, ...
        'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
    hold(ax_hdl, 'on');
end
camlight(ax_hdl)
ax_hdl.DataAspectRatio = [1,1,1];
view(ax_hdl, [1,1,1]);
ax_hdl.XLim = [0, size(labeled_array, 2)];
ax_hdl.YLim = [0, size(labeled_array, 1)];
ax_hdl.ZLim = [0, size(labeled_array, 3)];
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.Title.String = sprintf('#CC: %d', num_vis_label);
if nargout > 1
    varargout{1} = fig_hdl;
end
end