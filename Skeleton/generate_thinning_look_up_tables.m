%% Generate all the possible 3x3x3 binary cube configurations
num_neighbor_voxel = 27;
num_config = 2 ^ num_neighbor_voxel;
bin_config = logical(de2bi(0:num_config-1));
sk_bw = bwskel(image_mask);
% Remove all the configurations with center point equals false
bin_config = bin_config(bin_config(:,14)==true, :);
%% Construct look up table
% Remove end points
is_removable_point_Q = true(num_config/2,1);
% End point cannot be removed
is_endpoint_Q = (sum(bin_config,2)== 2);
is_removable_point_Q(is_endpoint_Q) = false;
% Do nothing if center point is false
% empty_center_point = (bin_config(:,14) == false);
% can_be_removed_table(empty_center_point) = false;
% Center point has all the 6 neighbor. 
is_internal_point_Q  = all(bin_config(:, [5, 11, 13, 15, 17, 23]),2);
is_removable_point_Q(is_internal_point_Q) = false;
% Find Euler-invariant points
is_Euler_invariant_removal_Q = fun_skeleton_is_Euler_invariant_removal(bin_config); 
% Find simple point
is_simple_point_removal_Q = fun_skeleton_is_simple_point_removal(bin_config);
% Modified version
is_removable_point_Q = is_removable_point_Q & is_Euler_invariant_removal_Q & is_simple_point_removal_Q;
save('./Skeleton/Skeleton_metadata/is_Euler_invariant_removal_Q.mat', 'is_Euler_invariant_removal_Q');
save('./Skeleton/Skeleton_metadata/is_simple_point_removal_Q.mat', 'is_simple_point_removal_Q');
save('./Skeleton/Skeleton_metadata/is_removable_point_Q.mat', 'is_removable_point_Q');