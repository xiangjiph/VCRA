clc;clear;
DataManager = FileManager;
dataset_name = 'WholeBrain';
image_grid_version = '240_cube';
skel_version = '240_cube_rec';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
allen_atlas = load('Allen_atlas.mat');
registration_version = 'Allen_2017_25um_landmark.mat';
Allen_atlas_id = load('Allen_atlas_id.mat');

pO2_folder_name = 'pO2_min_r_2um_sc';
reconstruction_version = '240_cube_recon_sc';
wb_stat_folder_name = 'whole_brain_stat_sc';
% Derivative variable
num_stack = numel(stack_list);
plot_stack_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
%%  Load registration and brain mask 
tic_load_data = tic;
wb_data_cell = cell(num_stack, 1);
for iter_stack = 1 : num_stack
    tmp_tic = tic;
    tmp_stack = stack_list{iter_stack};
    wb_data_cell{iter_stack} = fun_analysis_load_whole_brain_data_for_regional_analysis(...
        dataset_name, tmp_stack, image_grid_version, reconstruction_version, ...
        skel_version, registration_version, wb_stat_folder_name, true);
    
%     tmp_pO2_stat_filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, tmp_stack), ...
%         sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, tmp_stack, skel_version, pO2_folder_name));
%     wb_data_cell{iter_stack}.cube_stat.pO2_data = DataManager.load_data(tmp_pO2_stat_filepath);
    fprintf('Finish loading grid, registration and brain mask distance transform for stack %s. Elapsed time is %f seconds\n', tmp_stack, toc(tmp_tic));
end
fprintf('Finish loading data. Elapsed time is %f seconds.\n', toc(tic_load_data));
%% data extraction setting 
allen_atlas_map_old_id_to_new_id = wb_data_cell{1}.registration.map_oldID_to_newID;

de_opt_str = struct;
de_opt_str.cube_stat_Q = true;
de_opt_str.node_feature_Q = true;
de_opt_str.link_feature_Q = true;
de_opt_str.vessel_graph_Q = false;
de_opt_str.depth_dependent_density_Q = false;
de_opt_str.merge_cc_Q = false;
de_opt_str.save_Q = false;
de_opt_str.save_fp = [];
merge_stack_name = 'all_stack';
if de_opt_str.merge_cc_Q
    num_cc = 1;
else
    num_cc = 2;
end
%%
parent_str_name = 'Major_brain_regions';
analysis_region_id_list = {453 961 1080 477 549 [302 294] 4 771 354 528 776, 1097}; % Remove hypothalamus, use hippocampus region
vis_subfolder_name = 'Region_comparsion';

tmp_vis_region_list_ind = 1 : numel(analysis_region_id_list);

num_region = numel(analysis_region_id_list);
analysis_region_name = cell(num_region, 1);
for iter_region = 1 : num_region
    tmp_atlas_id = analysis_region_id_list{iter_region};
    if isscalar(tmp_atlas_id)
        analysis_region_name{iter_region} = allen_atlas.structure_table.name{full(allen_atlas.id_2_ind(tmp_atlas_id))};
    end
    tmp_atlas_id = full(allen_atlas_map_old_id_to_new_id(tmp_atlas_id));
    analysis_region_id_list{iter_region} = tmp_atlas_id;
end
analysis_region_name{6} = 'Superior colliculus';

visualization_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    merge_stack_name), vis_subfolder_name, parent_str_name);
%% Collect all the region statistics in all the stack
region_data = cell(num_region, num_stack);
collect_data_tic = tic;
for iter_stack = 1 : num_stack
    for iter_region = 1 : num_region
        tmp_tic = tic;
        tmp_region_id = analysis_region_id_list{iter_region};
        tmp_region_name = strrep(analysis_region_name{iter_region}, ' ', '_');
        
        region_data{iter_region, iter_stack} = fun_analysis_get_atlas_regional_data(wb_data_cell{iter_stack}, ...
            tmp_region_id, tmp_region_name, de_opt_str);
        toc(tmp_tic);
    end
end
fprintf('Finish collecting all the region data. Elapsed time is %f seconds.\n', toc(collect_data_tic));