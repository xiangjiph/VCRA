function feat_stat = fun_analysis_get_graph_feature_basic_stat(feature_table, validQ)
% fun_analysis_get_graph_feature_basic_stat compute the basic statistical
% properties for the features in the feature talbe that are numerical
% vector. 
% Input: 
%   feature_table: N-by-M table, where N is the number of object and M is
%   the number of feature components
%   validQ: N-by-1 logical vector. If provided, only components whose value
%   is true are used for computing the basic statistics. 
% Output: 
%   feat_stat: structure with fields, whose name are the computed feature
%   name in the given feature table. Each field is a structure returned by
%   fun_analysis_get_basic_statistics
% Note: 
% 1. Only the values in feature vector that are finite are used for
% computation. In other words, inf and nan are removed. 
%
if nargin < 2
    validQ = [];
end
if isempty(validQ)
    need_to_mask_Q = false;
else
    need_to_mask_Q = true;
end
feat_stat = struct;
if isa(feature_table, 'struct')
    feature_name_list = fieldnames(feature_table);
elseif isa(feature_table, 'table')
    feature_name_list = feature_table.Properties.VariableNames;
end
num_node_features = numel(feature_name_list);
one_item_tableQ = (size(feature_table, 1) == 1);
for iter_feature = 1 : num_node_features 
    tmp_feature_name = feature_name_list{iter_feature};
    tmp_data = feature_table.(tmp_feature_name);
    if islogical(tmp_data) % Convert the logical array to numerical array. 
        tmp_data = single(tmp_data);
    end
    if isnumeric(tmp_data) && isvector(tmp_data)
        if one_item_tableQ && ~isscalar(tmp_data) 
            % Skip the vector in the single-item table
            continue;
        end
        if need_to_mask_Q
            feat_stat.(tmp_feature_name) = fun_analysis_get_basic_statistics(tmp_data(validQ));
        else
            feat_stat.(tmp_feature_name) = fun_analysis_get_basic_statistics(tmp_data);
        end
    end
end
end