function [ax_hdl, varargout]= fun_vis_errorbar_shaded(x, y, ye, ax_hdl)

if nargin < 4
    ax_hdl = axes(figure);
else
    hold(ax_hdl, 'on');
end

assert(isvector(x) && isvector(y));
% Need to sort the value
% [y, y_ind] = sort(y, 'ascend');
% x = x(y_ind);
% ye = ye(y_ind);
x_vex = cat(2, x, flip(x));
y_vex = cat(2, y - ye, flip(y + ye));
is_valid_vex_Q = isfinite(x_vex) & isfinite(y_vex);
x_vex = x_vex(is_valid_vex_Q);
y_vex = y_vex(is_valid_vex_Q);
% y_vex = max(0, y_vex);

plt_hdl = plot(ax_hdl, x, y, 'LineWidth', 1.5);
hold(ax_hdl, 'on');
patch_hdl = fill(ax_hdl, x_vex(:), y_vex(:), plt_hdl.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% patch_hdl.FaceColor = plt_hdl.Color;
if nargout > 1
    varargout{1} = plt_hdl;
    varargout{2} = patch_hdl;
end
end