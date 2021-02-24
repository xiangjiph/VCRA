function [linear_fit_obj, varargout] = fun_vis_linear_fit_data_in_loglog(fit_data_x, fit_data_y, ax_hdl)

if nargin < 3
    fig_hdl = figure;
    ax_hdl = axes(fig_hdl);
else
    hold(ax_hdl, 'on');
    fig_hdl = ax_hdl.Parent;
end
% scatter(ax_hdl, fit_data_x, fit_data_y, 300, '.');
ax_hdl.XScale = 'log';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
valid_data_pair_Q = (fit_data_x > 0) & (fit_data_y > 0);
fit_data_x = fit_data_x(valid_data_pair_Q);
fit_data_y = fit_data_y(valid_data_pair_Q);
fit_data_x = log10(fit_data_x);
fit_data_y = log10(fit_data_y);
linear_fit_obj = fitlm(fit_data_x, fit_data_y);

hold(ax_hdl, 'on');
box(ax_hdl, 'on');
plt_hdl = plot(ax_hdl, 10 .^ fit_data_x, 10 .^ (linear_fit_obj.Coefficients.Estimate(1) + ...
    linear_fit_obj.Coefficients.Estimate(2) .* fit_data_x), 'LineWidth', 2);

ldg_hdl = legend(plt_hdl, sprintf('Exponent: %.2f \\pm %.2f\nIntercept: %.2f \\pm %.2f\nR-squared: %.3f\nData size: %d', ...
    linear_fit_obj.Coefficients.Estimate(2), linear_fit_obj.Coefficients.SE(2), ...
    linear_fit_obj.Coefficients.Estimate(1), linear_fit_obj.Coefficients.SE(1), ...
    linear_fit_obj.Rsquared.Adjusted, linear_fit_obj.NumObservations));

if nargout > 1
    fig_hdl.UserData.x = fit_data_x;
    fig_hdl.UserData.y = fit_data_y;
    varargout{1} = fig_hdl;
end
end