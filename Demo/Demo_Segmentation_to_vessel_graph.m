clc;clear;close all;
fprintf('Load demo data\n');
demo_data = load('Graph_refinemnt_test_data.mat');
vessel_image = demo_data.vessel_image;
vessel_mask = demo_data.vessel_mask;
%% Convert vessel mask to vessel graph: 
fprintf('Convert vessel mask to skeleton\n');
vessel_skeleton = bwskel(vessel_mask);
% Refine skeleton voxel position accoridng to image intensity
fprintf('Vessel skeleton centering\n');
vessel_skeleton_rc = fun_skeleton_recentering_within_mask(vessel_skeleton, vessel_image, vessel_mask);
% Convert vessel skeleton to graph
fprintf('Convert vessel skeleton to graph\n');
vessel_graph = fun_skeleton_to_graph(vessel_skeleton_rc);
% Add radius to the graph
fprintf('Compute distance transform with respect to the vessel mask to estimate the vessel radius\n');
vessel_mask_dt = bwdist(~vessel_mask);
vessel_graph = fun_graph_add_radius(vessel_graph, vessel_mask_dt);
fprintf('The vessel graph extracted from the vessel segmentation is shown below\n');
disp(vessel_graph)
% Reconstruct the vessel mask: 
fprintf('Generating 3D reconstruction of the vascular network from the graph\n');
vessel_mask_recon = fun_graph_reconstruction(vessel_graph);
volumeViewer(vessel_mask_recon)
