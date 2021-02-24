function [ax_hdl, varargout]= fun_vis_confidence_interval_shaded(x, y_med, yCI_l, yCI_h, ax_hdl)

if nargin < 5
    ax_hdl = axes(figure);
else
    hold(ax_hdl, 'on');
end

assert(isvector(x) && isvector(y_med));
if isrow(x)
    x = x.';
end
if isrow(y_med)
    y_med = y_med.';
end
if isrow(yCI_l)
    yCI_l = yCI_l.';
end
if isrow(yCI_h)
    yCI_h = yCI_h.';
end

is_valid_vex_Q = isfinite(x) & isfinite(y_med) & isfinite(yCI_l) & ...
    isfinite(yCI_h);
x = x(is_valid_vex_Q);
y_med = y_med(is_valid_vex_Q);
yCI_h = yCI_h(is_valid_vex_Q);
yCI_l = yCI_l(is_valid_vex_Q);
assert(all(yCI_l <= yCI_h), 'Upper confidence level should have value no less than the lower confidence level');
% Need to sort the value
x_vex = cat(2, x, flip(x));
y_vex = cat(2, yCI_l, flip(yCI_h));

plt_hdl = plot(ax_hdl, x, y_med, 'LineWidth', 1.5);
hold(ax_hdl, 'on');
patch_hdl = fill(ax_hdl, x_vex(:), y_vex(:), plt_hdl.Color, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
% patch_hdl.FaceColor = plt_hdl.Color;
if nargout > 1
    varargout{1} = plt_hdl;
    varargout{2} = patch_hdl;
end
end