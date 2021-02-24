function [vessel_graph, stat_str] = fun_analysis_whole_brain_graph_statistics(vessel_graph, ...
    dataset_name, stack, capillary_max_radius_um, capillary_min_radius_um, save_folder_name)

%%
persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end
if nargin < 6
    save_folder_name = 'whole_brain_stat';
end
save_im_folder = fullfile(DataManager.fp_analysis_data_folder(dataset_name, stack), ...
    save_folder_name);

stat_str = struct;
stat_str.dataset_name = dataset_name;
stat_str.stack = stack;
stat_str.fp = fullfile(save_im_folder, sprintf('%s_%s_%s.xml', dataset_name, stack, save_folder_name));
%% Compute basic features
if ~isfield(vessel_graph.link, 'features') || ~all(isfield(vessel_graph.link.features, {'dt_median', 'length'}))
    fprintf('Computing the radius and length of all the vessel segments\n');
    link_length = zeros(vessel_graph.link.num_cc, 1);
    link_med_r = zeros(vessel_graph.link.num_cc, 1);
    link_mean_r = zeros(vessel_graph.link.num_cc, 1);
    link_r_std = zeros(vessel_graph.link.num_cc, 1);
    tic
    for iter_link = 1 : vessel_graph.link.num_cc
        link_length(iter_link) = fun_graph_sub_to_length(fun_ind2sub(vessel_graph.num.mask_size, vessel_graph.link.cc_ind{iter_link}));
        tmp_r = full(vessel_graph.radius(vessel_graph.link.cc_ind{iter_link}));
        tmp_r = tmp_r(tmp_r >= capillary_min_radius_um);
        link_med_r(iter_link) = median(tmp_r);
        link_mean_r(iter_link) = mean(tmp_r);
        link_r_std(iter_link) = std(tmp_r);
    end
    toc
    vessel_graph.link.features.dt_mean = link_mean_r;
    vessel_graph.link.features.dt_median = link_med_r;
    vessel_graph.link.features.length = link_length;
    vessel_graph.link.features.dt_std = link_r_std;
end
%% Whole brain volume
brain_mask = DataManager.load_brain_mask(dataset_name, stack, 'whole_brain_d16x_annotated');
brain_mask_cc = bwconncomp(brain_mask > 0);
brain_mask_prop = regionprops(brain_mask_cc, 'BoundingBox');

stat_str.brain_volume_mm3 = numel(brain_mask_cc.PixelIdxList{1}) * 16^3 / 1e9;
fprintf('Volume of the brain is %f mm^3\n', stat_str.brain_volume_mm3);

stat_str.sample_size_mm = brain_mask_prop.BoundingBox(4:6) .* 16 / 1e3;
fprintf('Bounding box size: (%f, %f, %f) mm\n', stat_str.sample_size_mm);
%% Basic statistics
% Record statisitcs
stat_str.num.link = vessel_graph.link.num_cc;
stat_str.num.node = vessel_graph.node.num_cc;
stat_str.num.voxel = vessel_graph.num.skeleton_voxel;

has_no_endpoint_Q = all(vessel_graph.link.connected_node_label, 2);
has_1_endpoint_Q = any(vessel_graph.link.connected_node_label, 2) & ~has_no_endpoint_Q;
stat_str.num.link_wo_ep = nnz(has_no_endpoint_Q);
stat_str.num.link_w_1_ep = nnz(has_1_endpoint_Q);
stat_str.num.link_w_2_ep = nnz(~has_no_endpoint_Q & ~has_1_endpoint_Q);
fprintf('Number of links without endpoint: %d\n', stat_str.num.link_wo_ep);
stat_str.link.total_length_m = sum(vessel_graph.link.features.length, 'omitnan') / 1e6;
fprintf('Total link length: %f m\n', stat_str.link.total_length_m);
stat_str.link.wo_ep_total_length_m = sum(vessel_graph.link.features.length(has_no_endpoint_Q), 'omitnan') / 1e6;
fprintf('Total link without endpoint length: %f m\n', stat_str.link.wo_ep_total_length_m);

% Volume, length, surface and ratio
is_capillary_Q = vessel_graph.link.features.dt_median <= capillary_max_radius_um;
stat_str.capillary.total_length_m = sum(vessel_graph.link.features.length(is_capillary_Q), 'omitnan') / 1e6;
stat_str.capillary.length_ratio = stat_str.capillary.total_length_m ./ stat_str.link.total_length_m;
stat_str.link.total_surf_area_mm2 = sum(vessel_graph.link.features.length .* vessel_graph.link.features.dt_median * 2 * pi, 'omitnan') / 1e6;
stat_str.capillary.total_surf_area_mm2 = sum(vessel_graph.link.features.length(is_capillary_Q) .* vessel_graph.link.features.dt_median(is_capillary_Q) * 2 * pi, 'omitnan') / 1e6;
stat_str.capillary.surface_area_ratio = stat_str.capillary.total_surf_area_mm2 ./ stat_str.link.total_surf_area_mm2;
stat_str.link.total_volume_mm3 = sum(vessel_graph.link.features.length .* vessel_graph.link.features.dt_median.^2 * pi, 'omitnan') / 1e9;
stat_str.capillary.total_volume_mm3 = sum(vessel_graph.link.features.length(is_capillary_Q) .* vessel_graph.link.features.dt_median(is_capillary_Q) .^ 2 * pi, 'omitnan') / 1e9;
stat_str.capillary.volume_ratio = stat_str.capillary.total_volume_mm3 ./ stat_str.link.total_volume_mm3;

stat_str.link.length_density_m_mm3 = stat_str.link.total_length_m / stat_str.brain_volume_mm3;
stat_str.link.volume_density = stat_str.link.total_volume_mm3 / stat_str.brain_volume_mm3;

stat_str.capillary.length_density_m_mm3 = stat_str.capillary.total_length_m / stat_str.brain_volume_mm3;
stat_str.capillary.volume_density = stat_str.capillary.total_volume_mm3 / stat_str.brain_volume_mm3;
stat_str.capillary.max_radius_um = capillary_max_radius_um;
%% Number of skeleton voxel per cc without endpoints
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
hist_hdl = histogram(ax_hdl, vessel_graph.link.num_voxel_per_cc(has_no_endpoint_Q), 0:3:1000, 'Normalization', 'pdf');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Number of voxels';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.Title.String = 'Vessel segments without endpoints';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
stat_str.dist.link_wo_ep_num_vxl = fun_analysis_get_basic_statistics(vessel_graph.link.num_voxel_per_cc(has_no_endpoint_Q), true);
leg_hdl = legend(ax_hdl, fun_analysis_basic_stat_str_to_string(stat_str.dist.link_wo_ep_num_vxl));

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_link_wo_ep_voxel_num_distribution.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Node degree (self-loop counted once)
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
hist_hdl = histogram(ax_hdl, vessel_graph.node.num_link, 3:6, 'Normalization', 'pdf');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Degree';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.Title.String = 'Node degree';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
stat_str.dist.node_degree = fun_analysis_get_basic_statistics(vessel_graph.node.num_link, true);

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_link_voxel_num_distribution.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Find the connected component in the graph
num_link = vessel_graph.link.num_cc;
link_cc_label_list = fun_analysis_get_link_graph_cc_label(vessel_graph);
link_cc_label_cell = fun_bin_data_to_idx_list(link_cc_label_list);
num_link_per_cc = cellfun(@numel, link_cc_label_cell);
[num_link_largest_cc, giant_cc_id] = max(num_link_per_cc);
fprintf('The lartest connected component has %d links, which contains %f%% of the links\n', num_link_largest_cc, num_link_largest_cc / num_link * 100);

stat_str.graph.num_cc = numel(link_cc_label_cell);
stat_str.graph.cc_num_link = num_link_per_cc;
stat_str.graph.largest_cc_num_link = num_link_largest_cc;
stat_str.graph.largest_cc_link_fraction = num_link_largest_cc / num_link;

% Find the number of nodes in the largest connected component
giant_cc_link_label = double(link_cc_label_cell{giant_cc_id});
giant_cc_node_label = unique(full(vessel_graph.link.connected_node_label(giant_cc_link_label)));
stat_str.graph.largest_cc_num_node = numel(giant_cc_node_label);
stat_str.graph.largest_cc_node_fraction = stat_str.graph.largest_cc_num_node / vessel_graph.node.num_cc;
stat_str.graph.largest_cc_total_length = sum(vessel_graph.link.features.length(giant_cc_link_label));
is_capillary_in_giant_cc_Q = vessel_graph.link.features.dt_median(giant_cc_link_label) <= capillary_max_radius_um;
stat_str.graph.largest_cc_total_capillary_length = sum(vessel_graph.link.features.length(giant_cc_link_label(is_capillary_in_giant_cc_Q)));
stat_str.graph.largest_cc_capillary_length_ratio = stat_str.graph.largest_cc_total_capillary_length / stat_str.graph.largest_cc_total_length;
%% Link length distribution
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
hist_hdl = histogram(ax_hdl, vessel_graph.link.features.length, 0:2.5:1000, 'Normalization', 'pdf');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Link length (\mum)';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_link_length_distribution.png', dataset_name, stack, save_folder_name));
stat_str.dist.link_length = fun_analysis_get_basic_statistics(vessel_graph.link.features.length, true);
leg_hdl = legend(fun_analysis_basic_stat_str_to_string(stat_str.dist.link_length ));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Link length distribution - CDF
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
hist_edge = [10 .^ (linspace(log10(2), log10(1000), 30))];
hist_hdl = histogram(ax_hdl, vessel_graph.link.features.length, hist_edge, 'Normalization', 'cdf');
ax_hdl.XScale = 'log';
ax_hdl.YScale = 'linear';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Link length (\mum)';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_link_length_distribution_CDF.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

stat_str.link.length_90_prctile = prctile(vessel_graph.link.features.length, 90);
stat_str.link.length_10_prctile = prctile(vessel_graph.link.features.length, 10);
%% Link without endpoint length distribution
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
histogram(ax_hdl, vessel_graph.link.features.length(has_no_endpoint_Q), 0:2.5:1000, 'Normalization', 'pdf');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Link length (\mum)';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
stat_str.dist.link_wo_ep_length = fun_analysis_get_basic_statistics(vessel_graph.link.features.length(has_no_endpoint_Q));
leg_hdl = legend(fun_analysis_basic_stat_str_to_string(stat_str.dist.link_wo_ep_length));
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_link_wo_ep_length_distribution.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Capillary length distribution 
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
histogram(ax_hdl, vessel_graph.link.features.length(is_capillary_Q), 0:2.5:1000, 'Normalization', 'pdf');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Link length (\mum)';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
stat_str.dist.link_cap_length = fun_analysis_get_basic_statistics(vessel_graph.link.features.length(is_capillary_Q), true);
leg_hdl = legend(fun_analysis_basic_stat_str_to_string(stat_str.dist.link_cap_length));

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_capillary_length_distribution.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Radius distribution
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
hist_edge = [1 : 0.25 : 10, 10 .^ (linspace(log10(10.5), log10(150), 30))];
hist_hdl = histogram(ax_hdl, vessel_graph.link.features.dt_median, hist_edge, 'Normalization', 'pdf');
ax_hdl.XScale = 'log';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Radius (\mum)';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
stat_str.dist.link_radius = fun_analysis_get_basic_statistics(vessel_graph.link.features.dt_median, true);
leg_hdl = legend(fun_analysis_basic_stat_str_to_string(stat_str.dist.link_radius));
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_radius_distribution.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Radius distribution with fitting at multiple scale
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1.5];
ax_hdl = axes(fig_hdl);
hist_edge = [1 : 0.25 : 10, 10 .^ (linspace(log10(10.5), log10(100), 30))];
hist_hdl = histogram(ax_hdl, vessel_graph.link.features.dt_median, hist_edge, 'Normalization', 'pdf');
ax_hdl.XScale = 'log';
ax_hdl.YScale = 'log';
grid(ax_hdl, 'on');
ax_hdl.XLabel.String = 'Vessel segment radius (\mum)';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
% Fitting data between 3.5 and 50
[large_vessel_count, large_vessel_bin] = histcounts(vessel_graph.link.features.dt_median, 10 .^ (linspace(log10(capillary_max_radius_um), log10(50), 30)), 'Normalization', 'pdf');
large_vessel_r_val = movmean(large_vessel_bin, 2, 'Endpoints', 'discard');
[linear_fit_obj, fig_hdl] = fun_vis_linear_fit_data_in_loglog(large_vessel_r_val, ...
    large_vessel_count, ax_hdl);

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_seg_radius_distribution_with_noncap_scaling_fit.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Radius coefficient of variance 
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
dt_cv = vessel_graph.link.features.dt_std ./ vessel_graph.link.features.dt_mean;
hist_edge = 0 : 0.05 : 1;
histogram(ax_hdl, dt_cv, hist_edge, 'Normalization', 'cdf');
ax_hdl.XLabel.String = '\sigma_r / \mu_r';
ax_hdl.YLabel.String = 'CDF';
ax_hdl.FontSize = 14;
stat_str.dist.dt_std2mean = fun_analysis_get_basic_statistics(dt_cv, true);
leg_hdl = legend(ax_hdl, fun_analysis_basic_stat_str_to_string(stat_str.dist.dt_std2mean));
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_seg_radius_dt_std2mean.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Radius coefficient of variance 
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
dt_std2median = vessel_graph.link.features.dt_std ./ vessel_graph.link.features.dt_median;
hist_edge = 0 : 0.05 : 1;
histogram(ax_hdl, dt_std2median, hist_edge, 'Normalization', 'cdf');
ax_hdl.XLabel.String = '\sigma_r / median r';
ax_hdl.YLabel.String = 'CDF';
ax_hdl.FontSize = 14;
stat_str.dist.dt_std2median = fun_analysis_get_basic_statistics(dt_std2median, true);
leg_hdl = legend(ax_hdl, fun_analysis_basic_stat_str_to_string(stat_str.dist.dt_std2median));
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_seg_radius_dt_std2median.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Capillary radius coefficient of variance
fig_hdl = figure;
ax_hdl = axes(fig_hdl);
dt_std2median_cap = dt_std2median(is_capillary_Q);
hist_edge = 0 : 0.05 : 1;
histogram(ax_hdl, dt_std2median_cap, hist_edge, 'Normalization', 'cdf');
ax_hdl.XLabel.String = 'Capillary \sigma_r / median r';
ax_hdl.YLabel.String = 'CDF';
ax_hdl.FontSize = 14;
stat_str.dist.cap_dt_std2median = fun_analysis_get_basic_statistics(dt_std2median_cap, true);
leg_hdl = legend(ax_hdl, fun_analysis_basic_stat_str_to_string(stat_str.dist.cap_dt_std2median));
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_cap_radius_dt_std2median.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Link length - radius joint distribution
% Joint distribution
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, vessel_graph.link.features.length, vessel_graph.link.features.dt_median, 0 : 2.5 : 400,  10 .^ (0 : 0.1 : 2.5), 'DisplayStyle', 'tile', 'Normalization', 'pdf');
grid(ax_hdl, 'on');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'log';
ax_hdl.XLabel.String = 'Link length (\mum)';
ax_hdl.YLabel.String = 'Radius (\mum)';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';

ax_hdl.ColorScale = 'log';
cb = colorbar;
cb.Label.String = 'PDF';
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_link_length_radius_joint_dist.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Braching order distribution from large vessels
link_type = ones(vessel_graph.link.num_cc, 1);
link_type(~is_capillary_Q) = 2;
tic
link_order = fun_analysis_get_capillary_order_to_labeled_vessels(vessel_graph, link_type, 1, 2);
toc

valid_Q = ~isnan(link_order);
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
histogram(ax_hdl, link_order(valid_Q), 0.5:1:12.5, 'Normalization', 'pdf');
grid(ax_hdl, 'on');
ax_hdl.XScale = 'linear';
ax_hdl.YScale = 'linear';
ax_hdl.XLabel.String = 'Capillary branching order';
ax_hdl.YLabel.String = 'PDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
stat_str.dist.cap_branching_order = fun_analysis_get_basic_statistics(link_order(valid_Q), true);
leg_hdl = legend(fun_analysis_basic_stat_str_to_string(stat_str.dist.cap_branching_order));

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_%s_link_r_le_3500nm_branching_order.png', dataset_name, stack, save_folder_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
vessel_graph.link.features.capillary_branching_order = link_order;
end