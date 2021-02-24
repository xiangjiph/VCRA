function dt_mask = fun_graph_get_endpoint_link_mask(ep_ind, link_mask, valid_r)
% fun_graph_get_endpoint_link_mask generate a logical array where each
% point with value ture is closer to the endpoint than to any other link voxels
% in the link. 
% Input: 
%   ep_ind: numerical scaler, linear index of the position of the endpoint
%   link_mask: 3D logical array, mask of the link that the endpoing
%   belong to ( assume endpoint to be a member of the link voxel set) 
%   valid_r: numerical scalar. If the distance between a link voxel and the
%   endpoint is less than this value, they are exluced before the distance
%   transform. 
% Output: 
%   dt_mask: 3D logical array of the same size as link_mask. 
% Note: 
% 1. Here, the distance is the eculidean distance (default in bwdist). 
if nargin < 3
    valid_r = 3;
end
ep_mask = false(size(link_mask));
ep_mask(ep_ind) = true;
ep_dt = bwdist(ep_mask);

link_wo_ep_mask = link_mask & (ep_dt >= valid_r);
link_wo_ep_dt = bwdist(link_wo_ep_mask);


% Take the equal sign for rubusness issue. This function is mainly a
% subfunction for finding the linker between broken vessel segments. The
% boundary voxel might be equal-distant to more than 1 target link voxels.
% Using equal size can avoid the case where the boudnary voxel is excluded
% from the searching mask. 
dt_mask = ( ep_dt <= link_wo_ep_dt) ;
end