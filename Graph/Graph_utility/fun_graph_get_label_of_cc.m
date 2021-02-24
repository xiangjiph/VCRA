function cc_label_str = fun_graph_get_label_of_cc(map_ind_2_label, cc_ind_list, verboseQ)
% fun_graph_get_label_of_cc get the label of the connected component in the
% graph. 
% Input: 
%   map_ind_2_label: sparse matrix, field of vessel_graph.node/link
%   cc_ind_list: N-by-1 cell array of connected component indices in the mask 
% Output: 
%   cc_label: N-by-1 numerical vector, label of the connected component in
%   the graph 
if nargin < 3
    verboseQ = true;
end

num_kept_link = numel(cc_ind_list);
cc_label = cell(num_kept_link, 1);
contain_zero_Q = false(num_kept_link, 1);
num_label_per_cc = zeros(num_kept_link, 1);
for iter_link = 1 : num_kept_link
    tmp_label = full(map_ind_2_label(cc_ind_list{iter_link}));
    tmp_label_unique = unique(tmp_label);
    tmp_is_0_Q = (tmp_label_unique == 0);
    contain_zero_Q(iter_link) = any(tmp_is_0_Q);
    num_label_per_cc(iter_link) = numel(tmp_label_unique);
    cc_label{iter_link} = tmp_label_unique;
end
contain_more_than_1_label_Q = num_label_per_cc > 1;
unique_Q = ~contain_more_than_1_label_Q & ~contain_zero_Q;


label_list = cat(1, cc_label{:});
label_list = label_list(label_list > 0);

if ~all(num_label_per_cc == 1) && verboseQ
    warning('Link voxel in the graph correspond to more than one link in the new graph');
end
if any(contain_zero_Q) && verboseQ
    warning('Contains voxel not in the map')
end
cc_label_str.cc_label_cell = cc_label;

cc_label_str.unique_cc_label_list = cat(1, cc_label{unique_Q});
cc_label_str.cc_label_list = label_list;
cc_label_str.num_label_per_cc = num_label_per_cc;
cc_label_str.contain_zero_Q = contain_zero_Q;
cc_label_str.contain_more_than_1_label_Q = contain_more_than_1_label_Q;
end