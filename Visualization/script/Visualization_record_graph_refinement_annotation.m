DataManager = FileManager;
dataset_name = opt.grid_c.dataset_name;
stack = opt.grid_c.stack;
%%
annotation_group_name = 'link_ep1';
exp_label = 0;
im_save_folder = fullfile(DataManager.fp_visualization_folder(...
    dataset_name, stack), 'Graph_refinement', annotation_group_name);
%%
tmp_fig_hdl = gcf;
tmp_fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_%s_%d.png', ...
    dataset_name, stack, 'Graph_refinement', annotation_group_name, exp_label));
fun_print_image_in_several_formats(tmp_fig_hdl, tmp_fig_fp);
exp_label = exp_label + 1;
%% Dim short link 
annotation_group_name = 'dim_short_link';
exp_label = 0;
im_save_folder = fullfile(DataManager.fp_visualization_folder(...
    dataset_name, stack), 'Graph_refinement', annotation_group_name);
    %%
tmp_fig_hdl = gcf;
tmp_fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_%s_%d.png', ...
    dataset_name, stack, 'Graph_refinement', annotation_group_name, exp_label));
fun_print_image_in_several_formats(tmp_fig_hdl, tmp_fig_fp);
exp_label = exp_label + 1;
%% Linker
annotation_group_name = 'linker';
exp_label = 0;
im_save_folder = fullfile(DataManager.fp_visualization_folder(...
    dataset_name, stack), 'Graph_refinement', annotation_group_name);
%%
tmp_fig_hdl = gcf;
tmp_fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_%s_%d.png', ...
    dataset_name, stack, 'Graph_refinement', annotation_group_name, exp_label));
fun_print_image_in_several_formats(tmp_fig_hdl, tmp_fig_fp);
exp_label = exp_label + 1;
%%
 [tmp1, tmp2, tmp3] = fun_graph_get_linker_for_int_link(vessel_image, vessel_mask_1, vessel_skl, int_link_1);