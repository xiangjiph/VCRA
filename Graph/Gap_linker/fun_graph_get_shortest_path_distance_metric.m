function metric = fun_graph_get_shortest_path_distance_metric(int_array, mask_array)
% fun_graph_get_shortest_path_distance_metric generates a distance metric
% array for shortest path searching 
% Input: 
%   int_array: 3D numerical array. Intensity of voxel 
%   mask_array: 3D logical array. True for valid voxels
% Output: 
%   metric: 3D float array. 
% To do list: 
% 1. Consider how to integrate both the distance transform and the local
% itnensity - currently some linkers are too curvy due to the intensity
% advantage. E.g. region near the large bright vessels 
% 2. Consider how to integrate the orientation of the link with single
% endpoints. 
% 3. Prevent back tracking? 

% Parameters for mouselight image
min_neighbor_dist_cutoff = 10000;
local_int_max_cutoff = 20000;

if ~isfloat(int_array)
    int_array = single(int_array);
end

local_int_max = max(int_array(:));
local_int_min = min(int_array(:));
local_int_max = min(local_int_max_cutoff, local_int_max);
min_neighbor_dist = max(min_neighbor_dist_cutoff, (local_int_max - local_int_min)/2);
metric = max(0, local_int_max - int_array) + min_neighbor_dist;
metric(~mask_array) = inf;

end
