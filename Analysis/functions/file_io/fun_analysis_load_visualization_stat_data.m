function region_stat_str = fun_analysis_load_visualization_stat_data(dataset_name, stack_list, cc_list, compare_feature_group, ...
    compare_feature_subgroup, compare_structure_id, compare_task_name)

if nargin < 7
    compare_task_name = [];
end

persistent allen_atlas DataManager
if isempty(allen_atlas)
    allen_atlas = load('Allen_atlas.mat');
end
if isempty(DataManager)
    DataManager = FileManager;
end
compare_structure_name = allen_atlas.structure_table.name(allen_atlas.id_2_ind(compare_structure_id));
compare_structure_abv = allen_atlas.structure_table.acronym(allen_atlas.id_2_ind(compare_structure_id));
num_compare_group = numel(compare_feature_subgroup);

num_compare_cc = numel(cc_list);
num_stack = numel(stack_list);
num_compare_region = numel(compare_structure_id);
region_data = cell(num_compare_cc, num_stack, num_compare_region);

region_stat_str = struct;
region_stat_str.info.stack_list = stack_list;
region_stat_str.info.cc_list = cc_list;
region_stat_str.info.compare_structure_id = compare_structure_id;
region_stat_str.info.compare_structure_name = compare_structure_name;
region_stat_str.info.compare_structure_abv = compare_structure_abv;
region_stat_str.info.compare_feature_subgroup = compare_feature_subgroup;

for iter_group = 1 : num_compare_group
    tmp_subgroup_name = compare_feature_subgroup{iter_group};
    for iter_region = 1 : num_compare_region
        tmp_region_name = strrep(compare_structure_name{iter_region}, ' ', '_');
        for iter_stack = 1 : num_stack
            tmp_stack = stack_list{iter_stack};
            vis_data_root_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, tmp_stack), ...
                'cross_region_comparison');
            for iter_cc = 1 : num_compare_cc
                tmp_cc_name = cc_list{iter_cc};
                tmp_info = struct;
                tmp_info.dataset_name = dataset_name;
                tmp_info.stack = tmp_stack;
                tmp_info.compare_feature_group = compare_feature_group;
                tmp_info.compare_feature_subgroup = tmp_subgroup_name;
                tmp_info.cc_name = tmp_cc_name;
                tmp_info.task_name = compare_task_name;
                tmp_info.vis_data_folder = fullfile(vis_data_root_folder, compare_feature_group, tmp_cc_name, tmp_subgroup_name);
                tmp_info.region_name = tmp_region_name;
                tmp_info.region_name_wo_ds = compare_structure_name{iter_region};
                tmp_info.region_name_abv = compare_structure_abv{iter_region};
                tmp_data_root_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, tmp_stack), ...
                    compare_feature_group, tmp_region_name, tmp_cc_name);
                tmp_fp =  fullfile(tmp_data_root_folder, sprintf('%s_%s_%s_%s_%s.mat', tmp_stack, ...
                    tmp_cc_name, compare_feature_group, tmp_subgroup_name, tmp_region_name));
                tmp_data = DataManager.load_data(tmp_fp);
                tmp_info.data_field_name = fieldnames(tmp_data);
                tmp_data.info = tmp_info;
                region_data{iter_cc, iter_stack, iter_region} = tmp_data;
            end
        end
    end
    region_stat_str.(tmp_subgroup_name) = region_data;
    region_stat_str.plot_setting.(tmp_subgroup_name) = fun_vis_get_plot_setting(tmp_subgroup_name);
end
end

