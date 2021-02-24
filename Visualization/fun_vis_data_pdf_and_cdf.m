function [fig_hdl, varargout]= fun_vis_data_pdf_and_cdf(hist_edge, hist_pdf, hist_cdf, visibleQ)

if nargin < 4
    visibleQ = false;
end

if visibleQ
    fig_hdl = figure('Visible', 'on');
else
    fig_hdl = figure('Visible', 'off');
end
fig_hdl.Position(3:4) = fig_hdl.Position(3:4).* 2;
ax_hdl = axes(fig_hdl);
yyaxis(ax_hdl, 'right');
histogram(ax_hdl, 'BinCounts', hist_cdf, 'BinEdges', hist_edge);
ax_hdl.YLabel.String = 'CDF';
yyaxis(ax_hdl, 'left');
histogram(ax_hdl, 'BinCounts', hist_pdf, 'BinEdges', hist_edge);
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
ax_hdl.Title.String = strrep(plot_info_cell{2, iter_vis}, '_', ' ');
grid(ax_hdl, 'on');
leg_hdl = legend(tmp_leg_str, 'Location', 'best');
if nargout > 1
    varargout{1} = tmp_stat;
end
end