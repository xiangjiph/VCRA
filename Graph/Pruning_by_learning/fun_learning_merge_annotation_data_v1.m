function merged_data = fun_learning_merge_annotation_data(data, field_name)
% fun_learning_merge_annotation_data merges the annotation data
% % Input:
%   data: cell array of structures. Each structure contains several fields
%   including linker (for connecting the internal endpoint to the nearby
%   skeleton), link_ep1( link with 1 internal endpoint) and fake_link (link
%   form by segmentation error). Each field might contains one or more set
%   of data, which come from multiple round of graph refinement during the
%   annotation. For annotation pipeline, see
%   Automatic_annotation_pipeline_test_20190308.m. In each round of
%   annotation, there are 4 fields:
%       label: structure with fields to_removeQ and not_sureQ. Both are
%       logical array indicating whether the annotator decided this object
%       should be removed or not. Not_sureQ means the annotator was not
%       sure. 
%       raw_data: data from which the feature/data is generated from
%       feature/data: table of features that can be used for
%       classifier training
%       feature_name: generative field, name of the table column
% Output: 
%   merge_data: structure with fields
%       data: concatenated feautre table
%       label: structure with fields: 
%           removed
%
error('To be fix');
if ~isa(data, 'cell')
    data = {data};
end

num_block = numel(data);
merged_data = struct;
merged_data.data = {};
merged_data.label.to_removeQ = {};
merged_data.label.not_sureQ = {};
num_round_data = 0;
for idx = 1 : num_block
    tmp_data = data{idx};
    tmp_round_name = fieldnames(tmp_data.(field_name));
    tmp_num_round = numel(tmp_round_name);
    for iter_round = 1 : tmp_num_round
        num_round_data = num_round_data + 1;
        tmp_round_data = tmp_data.(field_name).(tmp_round_name{iter_round});
        merged_data.data{num_round_data} = tmp_round_data.data;
        merged_data.label.to_removeQ{num_round_data} = tmp_round_data.label.to_removeQ;
        merged_data.label.not_sure{num_round_data} = tmp_round_data.label.not_sureQ;
    end
end
merged_data.data = cat(1, merged_data.data{:});
merged_data.label.to_removeQ = cat(1, merged_data.label.to_removeQ{:});
merged_data.label.not_sureQ = cat(1, merged_data.label.not_sureQ{:});
end