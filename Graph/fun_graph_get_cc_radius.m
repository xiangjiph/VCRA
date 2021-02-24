function voxel_rad = fun_graph_get_cc_radius(cc_ind_list, mask_dt, method)
% fun_graph_get_cc_radius get the radius vector for the input link
% connected components. 
% Input: 
%   cc_ind_list: cell array, each one contains the voxel linear indices for
%   one connected component
%   mask_dt: numerical array, can be either sparse or full. If full, it is
%   the distance transform of the vessel mask. If sparse, it is the radius
%   field of the vessel graph. 
%   method: string. If 'all', use the median to replace 0 in the radius
%   readout( it is possible after gaps linking). If 'median', use the
%   median value for all the voxels in the link. 
% Output:
%   voxel_rad: numerical vector, of the same length as cat(1, cc_ind_list).
%   Radius of all the voxels in cc_ind_list, in the same order. 
%
% Implemented by Xiang Ji on 02/23/2019
%
if nargin < 3
    method = 'all';
end

num_cc = numel(cc_ind_list);
voxel_rad = cell(num_cc, 1);

if strcmp(method, 'median')
    use_median_Q = true;
elseif strcmp(method, 'all')
    use_median_Q = false;
end

if issparse(mask_dt)
    dt_is_sparse_Q = true;
else
    dt_is_sparse_Q = false;
end

for iter_cc = 1 : num_cc
    tmp_ind = cc_ind_list{iter_cc};
    if dt_is_sparse_Q
        tmp_dt = full(mask_dt(tmp_ind));
    else
        tmp_dt = mask_dt(tmp_ind);
    end
    tmp_dt_nonzero_Q = tmp_dt>0;
    if any(tmp_dt_nonzero_Q)
        tmp_median = median(tmp_dt(tmp_dt_nonzero_Q));
        if use_median_Q
            voxel_rad{iter_cc} = repelem(tmp_median, numel(tmp_ind), 1);
        else
            tmp_dt(~tmp_dt_nonzero_Q) = tmp_median;
            voxel_rad{iter_cc} = tmp_dt;
        end
    else
        voxel_rad{iter_cc} = repelem(0, numel(tmp_ind), 1);
    end
end
voxel_rad = cat(1, voxel_rad{:});
end