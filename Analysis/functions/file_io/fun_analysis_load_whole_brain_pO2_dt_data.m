function wb_pO2_data_str = fun_analysis_load_whole_brain_pO2_dt_data(grid_info, data_folder_name, layer_list)

%% Initialization
persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
dataset_name = grid_info.dataset_name;
stack = grid_info.stack;
grid_version = grid_info.version;
if nargin < 3
    layer_list = grid_info.layer;
end
num_core = 8;
%% Simplify the parfor loading
load_bbox_grid_label = grid_info.bbox_grid_label_array(:, :, layer_list);
load_bbox_grid_label = load_bbox_grid_label(load_bbox_grid_label > 0);
load_bbox_grid_ind = cat(1, grid_info.bbox_grid_ind{layer_list});
load_bbox_grid_sub = cat(1, grid_info.bbox_grid_sub{layer_list});
num_load_bbox = size(load_bbox_grid_label, 1);

[dt_lm_stat_cell, pO2_lm_stat_cell, local_dt_stat_cell, local_pO2_stat_cell, ...
    pO2_stat_in_dt_bin_cell, paired_extrema_cell] = deal(cell(num_load_bbox, 1));
% fit_krogh_coeff = nan(num_load_bbox, 1);
has_no_data_Q = false(num_load_bbox, 1);
%%
pool_obj = gcp('nocreate');
delete(pool_obj);
parpool(num_core);
parfor iter_bbox = 1 : num_load_bbox
	tmp_grid_label = load_bbox_grid_label(iter_bbox);
    % Load data
    try
        tmp_cube_data = DataManager.load_analysis_data_in_grid(data_folder_name, ...
            dataset_name, stack, grid_version, tmp_grid_label);
    catch ME
        has_no_data_Q(iter_bbox) = true;
        fprintf('Unable to load reconstructed mask (label: %d)\n', tmp_grid_label);
        fprintf('Error message: %s\n', ME.message);
        continue;
    end
    try
        % Precompute dt_lm and pO2_lm statistics in the simulation.         
        dt_lm_stat_cell{iter_bbox} = add_prctile(tmp_cube_data.dt_lm_stat, ...
            tmp_cube_data.local_dt_stat.val2ptrl_itp, ...
            tmp_cube_data.local_pO2_stat.val2ptrl_itp);
        
        pO2_lm_stat_cell{iter_bbox} = add_prctile(tmp_cube_data.pO2_lm_stat, ...
            tmp_cube_data.local_dt_stat.val2ptrl_itp, ...
            tmp_cube_data.local_pO2_stat.val2ptrl_itp);
        
        local_dt_stat_cell{iter_bbox} = tmp_cube_data.local_dt_stat;
        local_pO2_stat_cell{iter_bbox} = tmp_cube_data.local_pO2_stat;
        
        pO2_stat_in_dt_bin_cell{iter_bbox} = tmp_cube_data.pO2_stat_in_dt_bin;
        
%         fit_krogh_coeff(iter_bbox) = tmp_cube_data.fit_Krogh.corr_coeff;
        
%         tmp_cube_data.paired_extrema.pO2 = tmp_cube_data.pO2_lm.v(tmp_cube_data.paired_extrema.pO2_list_idx);
%         tmp_cube_data.paired_extrema.dt = tmp_cube_data.dt_lm.v(tmp_cube_data.paired_extrema.dt_list_idx);
%         paired_extrema_cell{iter_bbox} = rmfield(tmp_cube_data.paired_extrema, {'pO2_list_idx', 'dt_list_idx'});
    catch ME
        has_no_data_Q(iter_bbox) = true;
        fprintf('Unable to process cube (label: %d)\n', tmp_grid_label);
        rethrow(ME);
    end
end
pool_obj = gcp('nocreate');
delete(pool_obj);
%% Save to structure
wb_pO2_data_str = struct;
wb_pO2_data_str.grid_size = grid_info.grid_size;
wb_pO2_data_str.grid_label = load_bbox_grid_label;
wb_pO2_data_str.grid_ind = load_bbox_grid_ind;
wb_pO2_data_str.grid_sub = load_bbox_grid_sub;

wb_pO2_data_str.has_no_data_Q = has_no_data_Q;

wb_pO2_data_str.dt_lm = dt_lm_stat_cell;
wb_pO2_data_str.pO2_lm = pO2_lm_stat_cell;

wb_pO2_data_str.local_dt_stat = local_dt_stat_cell;
wb_pO2_data_str.local_pO2_stat = local_pO2_stat_cell;
wb_pO2_data_str.pO2_stat_in_dt_bin = pO2_stat_in_dt_bin_cell;
% wb_pO2_data_str.paired_extrema = paired_extrema_cell;
% wb_pO2_data_str.fit_krogh_coeff = fit_krogh_coeff;
end

function extrema_str = add_prctile(extrema_str, dt_to_prctile, pO2_to_prctile)

extrema_field_name = fieldnames(extrema_str);
num_field = numel(extrema_field_name);
for iter_field = 1 : num_field
    tmp_fn = extrema_field_name{iter_field};
    if contains(tmp_fn, 'dt')
       tmp_prctile = dt_to_prctile(extrema_str.(tmp_fn));
    elseif contains(tmp_fn, 'pO2')
        tmp_prctile = pO2_to_prctile(extrema_str.(tmp_fn));
    else
        disp(extrema_str);
        error('extrema_str does not contain dt and pO2 field');
    end
    extrema_str.(sprintf('%s_prctile', tmp_fn)) = tmp_prctile;
end


end