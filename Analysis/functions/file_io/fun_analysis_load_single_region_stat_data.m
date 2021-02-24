function region_data = fun_analysis_load_single_region_stat_data(dataset_name, stack_list, region_atlas_id)

persistent registration_str DataManager;
if isempty(registration_str)
    registration_str = load('Allen_atlas.mat');
end
if isempty(DataManager)
    DataManager = FileManager;
end

if isnumeric(region_atlas_id)
    region_list_ind = full(registration_str.id_2_ind(region_atlas_id));
    if region_list_ind > 0
        region_name = registration_str.structure_table.name(region_list_ind);
        region_name = region_name{1};
        region_name_abv = registration_str.structure_table.acronym(region_list_ind);
        region_name_abv = region_name_abv{1};
    else
       error('The structure atlas id does not exist'); 
    end
else
    region_name = region_atlas_id;
    region_name_abv = region_atlas_id;    
end
num_stack = numel(stack_list);
stack_data = cell(num_stack, 1);
valid_data_Q = false(num_stack, 1);
for iter_stack = 1 : num_stack
    tmp_stack = stack_list{iter_stack};
    region_data_folder = fullfile(DataManager.fp_metadata_folder(dataset_name, ...
        tmp_stack), 'region_data');
    region_data_filepath = fullfile(region_data_folder, ...
        sprintf('%s_%s_region_data_%s.mat', dataset_name, tmp_stack, ...
        strrep(region_name, ' ', '_')));
    if isfile(region_data_filepath)
        tmp_tic = tic;
        stack_data{iter_stack} = load(region_data_filepath);
        fprintf('Finish loading %s. Elapsed time is %f seconds\n', ...
            region_data_filepath, toc(tmp_tic));
        if isfield(stack_data{iter_stack}, 'cc_1') && isfield(stack_data{iter_stack}, 'cc_2')
            valid_data_Q(iter_stack) = true;
        end
    else
        fprintf('File %s does not exist\n', region_data_filepath);
    end
end

region_data = struct;
region_data.dataset_name = dataset_name;
region_data.stack_list = stack_list;
region_data.atlas_id = region_atlas_id;
region_data.region_name = region_name;
region_data.region_name_abv = region_name_abv;
region_data.data = stack_data;
region_data.valid_Q = valid_data_Q;
end