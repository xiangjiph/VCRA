function data_str = fun_vis_plot_multi_stack_histograms(stat_str_cell, opt)

if nargin < 2
    opt = [];
end
if isfield(opt, 'cdf_limit')
    plot_cdf_limit = opt.cdf_limit;
else
    plot_cdf_limit = [5e-3, 1 - 5e-3];
end
data_str = struct;

tmp_edge = cellfun(@(x) x.hist_edge, stat_str_cell, 'UniformOutput', false);
tmp_bin_val = cellfun(@(x) x.hist_bin_val, stat_str_cell, 'UniformOutput', false);
tmp_pdf = cellfun(@(x) x.hist_pdf, stat_str_cell, 'UniformOutput', false);

if isfield(opt, 'rebin_edge')
    for iter_cell = 1 : numel(tmp_edge)
        [tmp_pdf{iter_cell}, tmp_edge{iter_cell}] = ...
            fun_analysis_rebin_histcount_edge(tmp_pdf{iter_cell}, tmp_edge{iter_cell}, ...
            opt.rebin_edge);
        tmp_bin_val{iter_cell} = movmean(tmp_edge{iter_cell}, 2, 'Endpoints', 'discard');
    end
end

if isfield(opt, 'int_x')
    tmp_int = fun_analysis_get_xy_curve_avgNstd_by_interpolation(tmp_bin_val, tmp_pdf, opt.int_x);
else
    tmp_int = fun_analysis_get_xy_curve_avgNstd_by_interpolation(tmp_bin_val, tmp_pdf);
end

if isfield(opt, 'Figure')
    fig_hdl = opt.Figure;
else
    fig_hdl = figure;
end

if isfield(opt, 'Axes')
    ax_hdl = opt.Axes;
else
    ax_hdl = axes(fig_hdl);
end
if opt.HistQ
    for iter_cell = 1 : numel(stat_str_cell)
        [tmp_plot_count, tmp_plot_edge] = fun_analysis_select_histcount_edge_by_percentile(tmp_pdf{iter_cell}, ...
            tmp_edge{iter_cell}, plot_cdf_limit, 'pdf');
        scatter(ax_hdl, movmean(tmp_plot_edge, 2, 'Endpoints', 'discard'), ...
            tmp_plot_count, 30, 'filled', 'MarkerFaceAlpha', 1);
%         histogram(ax_hdl, 'BinEdges', tmp_plot_edge, 'BinCounts', tmp_plot_count);
        hold(ax_hdl, 'on');
    end
end
if opt.ErrorBarQ
    [ax_hdl, eb_line_hdl, eb_shade_hdl] = fun_vis_errorbar_shaded(tmp_int.interpolate_x, ...
        tmp_int.y_avg, tmp_int.y_std, ax_hdl);
    data_str.errorbar_hdl = eb_line_hdl;
    data_str.errorbar_patch_hdl = eb_shade_hdl;
%     errorbar(ax_hdl, tmp_int.interpolate_x, tmp_int.y_avg, tmp_int.y_std, 'LineWidth', 2);
end

if isfield(opt, 'XLim')
    ax_hdl.XLim = opt.XLim;
else
%     ax_hdl.XLim = tmp_plot_edge([1, end]);
end

if isfield(opt, 'YLim')
    ax_hdl.YLim = opt.YLim;
end

if isfield(opt, 'LegendString')
    leg_hdl = legend(ax_hdl, opt.LegendString);
    if isfield(opt, 'LegendLabel')
        leg_hdl.Label.String = opt.LegendLabel;
    end
    data_str.legend_hdl = leg_hdl;
end

if isfield(opt, 'TitleString')
    ax_hdl.Title.String = opt.TitleString;
end
if isfield(opt, 'XLabelString')
    ax_hdl.XLabel.String = opt.XLabelString;
end
if isfield(opt, 'YLabelString')
    ax_hdl.YLabel.String = opt.YLabelString;
end
if isfield(opt, 'YScale')
    ax_hdl.YScale = opt.YScale;
end
if isfield(opt, 'XScale')
    ax_hdl.XScale = opt.XScale;
end
%% Output

data_str.fig_hdl = fig_hdl;
data_str.ax_hdl = ax_hdl;
data_str.interpolation_str = tmp_int;

end