clc;clear;
DataManager = FileManager;
dataset_name = 'WholeBrain';
image_grid_version = '240_cube';
reconstruction_version = '240_cube_recon';
skel_version = '240_cube_rec';
stack_list = {'ML_2018_08_15', 'ML20190124'};
num_stack = numel(stack_list);
allen_atlas = load('Allen_atlas.mat');
registration_version = 'Allen_2017_25um_landmark.mat';
Allen_atlas_id = load('Allen_atlas_id.mat');

linear_scaling_factor = 1.0521;
%%  Load registration and brain mask 
wb_data_cell = cell(num_stack, 1);
tic_load_data = tic;
for iter_stack = 1 : num_stack
    wb_data_cell{iter_stack} = fun_analysis_load_whole_brain_data_for_regional_analysis(...
        dataset_name, stack_list{iter_stack}, image_grid_version, reconstruction_version, ...
        skel_version, registration_version);
end
fprintf('Finish loading data. Elapsed time is %f seconds.\n', toc(tic_load_data));
%% data extraction setting 
de_opt_str = struct;
de_opt_str.cube_stat_Q = true;
de_opt_str.node_feature_Q = true;
de_opt_str.link_feature_Q = true;
de_opt_str.vessel_graph_Q = true;
de_opt_str.depth_dependent_density_Q = true;
de_opt_str.merge_cc_Q = false;
de_opt_str.save_Q = false;
de_opt_str.save_fp = [];
merge_stack_name = 'all_stack';

allen_atlas_map_old_id_to_new_id = wb_data_cell{1}.registration.map_oldID_to_newID;
%%
parent_str_name = 'Cerebral_cortex_transition';
analysis_region_id_list = [453 500 247 669 895 31 961];
tmp_vis_region_list_ind = 1 : numel(analysis_region_id_list);
% analysis_region_id_list = analysis_region_id_list(tmp_vis_region_list_ind);

num_region = numel(analysis_region_id_list);
analysis_region_name = allen_atlas.structure_table.name(full(allen_atlas.id_2_ind(analysis_region_id_list)));
analysis_region_name_abbrv = allen_atlas.structure_table.acronym(full(allen_atlas.id_2_ind(analysis_region_id_list)));
vis_subfolder_name = 'Depth_dependence';
visualization_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    merge_stack_name), vis_subfolder_name, parent_str_name);

analysis_region_id_list = full(allen_atlas_map_old_id_to_new_id(analysis_region_id_list));
%% Collect all the region statistics in all the stack
region_data = cell(num_region, num_stack);
collect_data_tic = tic;
for iter_stack = 1 : num_stack
    for iter_region = 1 : num_region
        tmp_tic = tic;
        tmp_region_id = analysis_region_id_list(iter_region);
        tmp_region_name = strrep(analysis_region_name{iter_region}, ' ', '_');
        
        region_data{iter_region, iter_stack} = fun_analysis_get_atlas_regional_data(wb_data_cell{iter_stack}, ...
            tmp_region_id, tmp_region_name, de_opt_str);
        toc(tmp_tic);
    end
end
fprintf('Finish collecting all the region data. Elapsed time is %f seconds.\n', toc(collect_data_tic));
%% Capillary length density for individual region 
vis_target = 'vessel';
switch vis_target
    case 'capillary'
        y_label = 'Capillary length density (m/mm^3)';
        vis_subfolder_name = 'capillary_length_density';
    case 'vessel'
        y_label = 'Vessel length density (m/mm^3)';
        vis_subfolder_name = 'vessel_length_density';
end

cap_density_vs_depth_cell = cell(num_region, 1);
cc_label = {'Left', 'Right'};
for iter_region = 1 : num_region
    str_data_cell = cell(2, num_stack);

    cap_den_depth_data = struct;
    cap_den_depth_data.length_density = cell(2, num_stack);
    cap_den_depth_data.depth = cell(2, num_stack);
    cap_den_depth_data.legend = cell(2, num_stack);
    cap_den_depth_data.data_size = numel(cap_den_depth_data.depth);
    for iter_stack = 1 : num_stack
        for iter_cc = 1 : 2
            tmp_cc_data = region_data{iter_region, iter_stack}.(sprintf('cc_%d', iter_cc));
            str_data_cell{iter_cc, iter_stack} = tmp_cc_data;
            switch vis_target
                case 'capillary'
                    cap_den_depth_data.length_density{iter_cc, iter_stack} = tmp_cc_data.dt.cap_length_density ./ linear_scaling_factor^2;
                case 'vessel'
                    cap_den_depth_data.length_density{iter_cc, iter_stack} = tmp_cc_data.dt.ves_length_density ./ linear_scaling_factor^2;
            end
            cap_den_depth_data.depth{iter_cc, iter_stack} = tmp_cc_data.dt.bin_val .* linear_scaling_factor;
            cap_den_depth_data.legend{iter_cc, iter_stack} = sprintf('%s %s', strrep(stack_list{iter_stack}, '_', ''), cc_label{iter_cc});
        end
    end
    cap_den_depth_data.avg_interpolation = fun_analysis_get_xy_curve_avgNstd_by_interpolation(...
        cap_den_depth_data.depth, cap_den_depth_data.length_density);
    % Plot
    fig_hdl = figure('Visible', 'on');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 1.5;
    ax_hdl = axes(fig_hdl); %#ok<LAXES>
    hold(ax_hdl, 'on');
    for iter_cell = 1 : cap_den_depth_data.data_size
        plot(ax_hdl, cap_den_depth_data.depth{iter_cell}, cap_den_depth_data.length_density{iter_cell}, 'LineWidth', 1, 'LineStyle', '--');
    end       
    errorbar(ax_hdl, cap_den_depth_data.avg_interpolation.interpolate_x, ...
        cap_den_depth_data.avg_interpolation.y_avg, ...
        cap_den_depth_data.avg_interpolation.y_std, 'LineWidth', 3);
    tmp_legend = cat(1, cap_den_depth_data.legend(:), {'Average \pm STD'});
    leg_hdl = legend(tmp_legend{:}, 'Location', 'best');
    ax_hdl.FontSize = 14;
    ax_hdl.FontWeight = 'bold';
    box(ax_hdl, 'off');
    grid(ax_hdl, 'off');
    ax_hdl.XLabel.String = 'Distance to the cortical surface(\mum)';
    ax_hdl.YLabel.String = y_label;
    ax_hdl.Title.String = sprintf('%s', analysis_region_name{iter_region});
    ax_hdl.YLim(1) = 0;
    fig_name = fullfile(visualization_folder, vis_subfolder_name, sprintf('%s_%s_%s_vs_depth_in_%s.png', ...
        dataset_name, merge_stack_name, vis_subfolder_name, strrep(analysis_region_name{iter_region}, ' ' , '_')));
    fun_print_image_in_several_formats(fig_hdl, fig_name);
    delete(fig_hdl);
    cap_density_vs_depth_cell{iter_region} = cap_den_depth_data;
end
%% Capillary length density - show selected regions together
fig_hdl = figure('Visible', 'on');
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 1.5;
ax_hdl = axes(fig_hdl);
hold(ax_hdl, 'on');
leg_target_hdl = [];
for iter_region = 1 : num_region
%     if any(iter_region == tmp_vis_region_list_ind)
        tmp_x = cap_density_vs_depth_cell{iter_region}.avg_interpolation.interpolate_x;
        tmp_selected_Q = (tmp_x <= 1e3);
        if iter_region == 7
            tmp_selected_Q = (tmp_x <= 6e2);
        end        
        tmp_y = cap_density_vs_depth_cell{iter_region}.avg_interpolation.y_avg;
        tmp_y_err = cap_density_vs_depth_cell{iter_region}.avg_interpolation.y_std;
        tmp_x = tmp_x(tmp_selected_Q);
        tmp_y = tmp_y(tmp_selected_Q);
        tmp_y_err = tmp_y_err(tmp_selected_Q);
        [ax_hdl, tmp_plt_hdl, tmp_patch_hdl] = fun_vis_errorbar_shaded(...
            tmp_x, tmp_y, tmp_y_err, ax_hdl);
        tmp_plt_hdl.LineWidth = 3;
        tmp_patch_hdl.FaceAlpha = 0.3;
        leg_target_hdl = [leg_target_hdl, tmp_plt_hdl];
        if iter_region == 7
            tmp_plt_hdl.Color = 'k';
            tmp_patch_hdl.FaceColor = 'k';
        end
%         errorbar(ax_hdl, tmp_x, tmp_y, tmp_y_err, 'LineWidth', 3);
%         plot(ax_hdl, tmp_x, tmp_y, 'LineWidth', 3);
%     end
end
ax_hdl.YLim(1) = 0;
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
ax_hdl.XLabel.String = 'Distance to the cortical surface (\mum)';
ax_hdl.YLabel.String = y_label;
leg_hdl = legend(leg_target_hdl, analysis_region_name, 'Location', 'best');
leg_hdl.Title.String = 'Mean \pm STD';
fig_name = fullfile(visualization_folder, vis_subfolder_name, ...
    sprintf('%s_%s_%s_vs_depth_avg_selected_from_%s.png', dataset_name, merge_stack_name, vis_subfolder_name, parent_str_name));
fun_print_image_in_several_formats(fig_hdl, fig_name);