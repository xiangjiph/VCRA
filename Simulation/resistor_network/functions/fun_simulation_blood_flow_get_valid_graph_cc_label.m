function [valid_cc_label, varargout] = fun_simulation_blood_flow_get_valid_graph_cc_label(graph_cc_label, cc_type)
% fun_simulation_blood_flow_get_valid_graph_cc_label finds the graph
% connected components that have both artery (label 2) and vein (3), which
% can be used for blood flow simulation later. Notice that if more than 1
% valid labels are found, the adjaceny matrix will be block diagonal and
% the simulation can be done separately. 
% Input: 
%   graph_cc_label: numerical vector, graph connected component label of
%   the vessel link. if graph_cc_label(i) = j, it means the i-th vessel
%   link is in the graph connected component j. 
%   cc_type: numerical vector, type of vessel link. By default, 1 for
%   capillary / unlabeled, 2 for artery and 3 for vein 
% 
% Output: 
%   valid_cc_label: label of the graph connected components that contain
%   both artery and vein
%   graph_link_cc: additional information
% 
% Implemented by Xiang Ji on 07/25/2019

%% Parameters
artery_label = 2;
vein_label = 3;
%% 
[graph_link_cc.idx, graph_link_cc.cc_label] = fun_bin_data_to_idx_list(graph_cc_label);
graph_link_cc.num_cc = numel(graph_link_cc.idx);
graph_link_cc.contains_artery_and_vein_Q = false(graph_link_cc.num_cc, 1);
for iter_cc = 1 : graph_link_cc.num_cc
    if ~isnan(graph_link_cc.cc_label(iter_cc))
        if any(cc_type(graph_link_cc.idx{iter_cc}) == artery_label) && ...
                any(cc_type(graph_link_cc.idx{iter_cc}) == vein_label)
            graph_link_cc.contains_artery_and_vein_Q(iter_cc) = true;
        end
    end
end
fprintf('Number of valid cc for simulation is %d\n', nnz(graph_link_cc.contains_artery_and_vein_Q));
valid_cc_label = graph_link_cc.cc_label(graph_link_cc.contains_artery_and_vein_Q);

if nargout > 1
    varargout{1} = graph_link_cc;
end
end