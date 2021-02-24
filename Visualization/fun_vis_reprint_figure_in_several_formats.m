function output_fig_fp = fun_vis_reprint_figure_in_several_formats(fig_fp)
fig_hdl = openfig(fig_fp);
output_fig_fp = strrep(fig_fp, '_data.fig', '.png');
fun_print_image_in_several_formats(fig_hdl, output_fig_fp);
end