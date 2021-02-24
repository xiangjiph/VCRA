function pcl_str = fun_analysis_percolation_bond_rm_simulation(graph_str, bond_rm_p, num_simulation)

% Parameters
num_save_cc = 10;
%
pcl_str = struct;
pcl_str.bond_rm_prob = bond_rm_p;
pcl_str.bond_ocp_prob = 1 - bond_rm_p;
pcl_str.num_simulation = num_simulation;

[pcl_str.num_cc, pcl_str.largest_cc_num_site, pcl_str.num_site, ...
    pcl_str.largest_cc_fraction] = deal(nan(num_simulation, 1));
[pcl_str.cc_fraction]= deal(nan(num_save_cc, num_simulation));

num_bond = graph_str.numedges;

pcl_str.num_bond_0 = num_bond;
pcl_str.num_site_0 = graph_str.numnodes;
pcl_str.num_site_per_cc_0 = cellfun(@numel, conncomp(graph_str, 'OutputForm', 'cell'));
pcl_str.num_cc_0 = numel(pcl_str.num_site_per_cc_0);

num_rm_link = round(num_bond * bond_rm_p);
for iter_simu = 1 : num_simulation
    %% Randomly remove links, as well as newly created endpoints
    tmp_rm_link_label = randsample(num_bond, num_rm_link, false);
    tmp_graph = rmedge(graph_str, tmp_rm_link_label);
    tmp_site_graph_label = conncomp(tmp_graph, 'OutputForm', 'cell');
    %% Record graph connected component data
    tmp_cc_num_site = sort(cellfun(@numel, tmp_site_graph_label), 'descend');
    pcl_str.num_cc(iter_simu) = numel(tmp_cc_num_site);
    tmp_ind = 1 : min(num_save_cc, pcl_str.num_cc(iter_simu));
    pcl_str.num_site(iter_simu) = sum(tmp_cc_num_site);
    pcl_str.cc_fraction(tmp_ind, iter_simu) = tmp_cc_num_site(tmp_ind) / pcl_str.num_site(iter_simu);
    
    pcl_str.largest_cc_num_site(iter_simu) = tmp_cc_num_site(1);
    pcl_str.largest_cc_fraction(iter_simu) = pcl_str.cc_fraction(1, iter_simu);
end
pcl_str.cc_fraction = pcl_str.cc_fraction';
end