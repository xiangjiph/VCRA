function TD = fun_learning_save_annotation_result(label_name, label, feature, raw_data, TD, overwrite_field, saveQ)
% 
if nargin < 5
    TD = [];
    saveQ = false;
    overwrite_field = false;
elseif nargin < 6
    saveQ = false;
    overwrite_field = false;
elseif nargin < 7
    saveQ = false;    
end

tmp_str = fun_learning_get_training_data_from_annotation_label(label, feature, raw_data, false);
tmp_str_field_names = fieldnames(tmp_str);
field_name = cell(size(tmp_str_field_names));
num_filed_name = numel(field_name);
if ~isfield(TD, label_name)
    TD.(label_name) = [];
end

for iter_field = 1 : num_filed_name
    % Create field if not exist. 
    if ~isfield(TD.(label_name), tmp_str_field_names{iter_field}) || overwrite_field
        TD.(label_name).(tmp_str_field_names{iter_field}) = [];
    end
    % Append the structure to the existing cell
    TD.(label_name).(tmp_str_field_names{iter_field}){end + 1} = tmp_str.(tmp_str_field_names{iter_field});    
end
% if ~isfield(TD.(label_name), 'all_label') || overwrite_field
%     TD.(label_name).all_label = [];
% end
% TD.(label_name).all_label{end + 1} = label;

if saveQ
    if isfield(TD, 'info') && isfield(TD.info, 'filepath')
        folder_path = fileparts(TD.info.filepath);
        if ~isfolder(folder_path)
            mkdir(folder_path)
        end
        save(TD.info.filepath, '-struct', 'TD');
    end
end
end