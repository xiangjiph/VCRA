function exit_code = fun_analysis_compute_histogram_of_region_stat(region_stat, merge_dim)

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end

if nargin < 2
    merge_dim = [];
end
region_stat.info.folder_name = strrep(region_stat.info.region_name, ' ', '_');

plot_field_name = fieldnames(region_stat.plot_setting);
num_plot_field = numel(plot_field_name);

num_merge_dim = numel(merge_dim);

for iter_plot_field = 1 : num_plot_field
    tmp_plot_field_name = plot_field_name{iter_plot_field};
    tmp_data = region_stat.(tmp_plot_field_name);
    
    tmp_plot_info = region_stat.plot_setting.(plot_field_name{iter_plot_field});
    
    assert(all(cellfun(@istable, tmp_data), 'all'), 'All the data inside the cell array should be table');
    
    if num_merge_dim > 0
        tmp_data = fun_concatenate_cell_array_in_arbitrary_dims(tmp_data, merge_dim, 1);
        for iter_merge = 1 : num_merge_dim
            region_stat.info.dim_name_list{merge_dim(iter_merge)} = {sprintf('all_%s', ...
                region_stat.info.dim_name{merge_dim(iter_merge)})};
        end
    end
    remain_data_size = size(tmp_data);
    remain_data_numel = numel(tmp_data);
    
    cc_name_list = region_stat.info.dim_name_list{1};
    stack_name_list = region_stat.info.dim_name_list{2};
    
    for iter_cell = 1 : remain_data_numel
        assert(ismatrix(tmp_data), 'To be implemented');
        % Assuming two fields need to fix later
        [dim1_idx, dim2_idx] = ind2sub(remain_data_size, iter_cell);
        tmp_stack_name = stack_name_list{dim2_idx};
        tmp_cc_name = cc_name_list{dim1_idx};
        
        root_folder = fullfile(DataManager.fp_visualization_folder(region_stat.info.dataset_name, ...
            tmp_stack_name), region_stat.info.str_name, region_stat.info.folder_name, ...
            tmp_cc_name);
        plot_data = tmp_data{iter_cell};
        if isempty(plot_data)
            fprintf('%s is an empty dataset. Skip\n', root_folder);
            break;
        end  
        tmp_name_string = sprintf('%s_%s_%s_%s_%s', tmp_stack_name, tmp_cc_name, ...
            region_stat.info.str_name, tmp_plot_field_name, region_stat.info.folder_name);
        tmp_analysis_field = numel(tmp_plot_info);
        stat_str = struct;
        stat_str.info.stack = tmp_stack_name;
        stat_str.info.cc = tmp_cc_name;
        stat_str.info.process_date = datestr(now);
        stat_str.info.filepath = fullfile(root_folder, sprintf('%s.xml', tmp_name_string));
        stat_str.info.structure_name = region_stat.info.region_name;
        for iter_plot = 1 : tmp_analysis_field
            tmp_plot_feature_name = tmp_plot_info{iter_plot};
            tmp_plot_data = plot_data.(tmp_plot_feature_name);
            tmp_stat = fun_analysis_get_basic_statistics(tmp_plot_data);
            tmp_stat.stack = tmp_stack_name;
            tmp_stat.cc_name = tmp_cc_name;
            tmp_stat.region_name =  region_stat.info.region_name;
            % Save statistical data
            stat_str.(tmp_plot_feature_name) = tmp_stat;
        end
        fprintf('Finish writing PDF of region statistics features\n');
        xml_write(stat_str.info.filepath, stat_str);
        stat_str.info.filepath = strrep(stat_str.info.filepath, '.xml', '.mat');
        save(stat_str.info.filepath, '-struct', 'stat_str');        
        fprintf('Finish writing region statistics XML file\n');
    end
end
exit_code = 0;
end
