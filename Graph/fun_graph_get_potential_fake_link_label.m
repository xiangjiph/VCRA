function check_label = fun_graph_get_potential_fake_link_label(input_graph)
% This is an ad hoc function for finding the link that are potentially
% form by segmentation/ stitching artefacts. These links can be check by
% App_check_fake_gap_linker
%
% Parameters:
large_vessel_r = 5;
th_max_dt_fake_link_between_large_vessels = 2;
th_ep_dt_sum_2_length = 0.3;
th_short_loop_length = 100;
th_max_num_edges_in_loop = 5;
th_max_length_fake_link_between_capillaries = 40;
%% Ad hoc selection
% Links in short loops
short_loop_Q = input_graph.link.features.shortest_loop_lenght < th_short_loop_length;
% Links in loops of less than 6 links 
short_geodesic_loop_Q = input_graph.link.features.shortest_loop_geodesic_lenght <= th_max_num_edges_in_loop;

% The sum of the distance transform of the two end of the link is greater
% than the length of the link 
large_ep2length_Q = (input_graph.link.features.dt_ep_sum_2_ep_dist) > th_ep_dt_sum_2_length;

% Fake link between two large vessels ( happen frequently on the surface) 
link_between_large_vessels_Q = (input_graph.link.features.dt_ep1 >= large_vessel_r & ...
    input_graph.link.features.dt_ep2 >= large_vessel_r) & ...
    input_graph.link.features.dt_min <= th_max_dt_fake_link_between_large_vessels & large_ep2length_Q;

% Fake link between large vessels and surrounding vessels
link_between_large_vessel_and_neighbors_Q = ...
    (input_graph.link.features.dt_ep1 >= large_vessel_r | ...
    input_graph.link.features.dt_ep2 >= large_vessel_r) & ...    
    input_graph.link.features.dt_min <= 2 & large_ep2length_Q;

% Short link in short loops - to find artefact connection form between
% large vessels. Happens in cerebellum
short_link_in_short_loop_Q = short_geodesic_loop_Q & ...
    input_graph.link.features.length < 10 & ...
    large_ep2length_Q & ...
    input_graph.link.features.dt_max >= 3;
% Dim short link between capillaries
dim_short_link_Q = input_graph.link.features.int_min < 22000 & ...
    large_ep2length_Q & input_graph.link.features.length < th_max_length_fake_link_between_capillaries;
% Stitching artefact: straight line on the xy plane
stitching_artifact_line_Q = (input_graph.link.features.dt_min <= sqrt(2.1)) & ...
    (input_graph.link.features.length < 50) & ...
    (input_graph.link.features.shortest_loop_lenght < 150) & ...
    input_graph.link.features.not_near_boundary_Q;
% 
link_in_loop_of_large_vessels_Q = input_graph.link.features.dt_median > 3 &...
    short_loop_Q;
% 
not_length_1_Q = input_graph.link.features.length > 1;

tmp_check_link_Q = (short_loop_Q | link_between_large_vessels_Q | link_between_large_vessel_and_neighbors_Q | ...
    short_link_in_short_loop_Q | dim_short_link_Q | stitching_artifact_line_Q | link_in_loop_of_large_vessels_Q);
check_label = find(tmp_check_link_Q & not_length_1_Q & input_graph.link.features.not_near_boundary_Q);
end