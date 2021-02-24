function dt_pur_str = fun_analysis_DT_single_link_perturbation(vessel_label_mask, ...
    vessel_graph, connectivity_graph, link_label_list)

if nargin < 4
    link_label_list = 1 : vessel_graph.link.num_cc;
end
num_test_label = numel(link_label_list);
assert(isfield(vessel_graph.link.features, {'nearest_dt_ind'}));
%% Initialization








end