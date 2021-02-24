function [ori_vec, varargout]= fun_analysis_get_ordered_vxl_sub_list_ori_vec(voxel_sub, method)

if nargin < 2
    method = 'ep2ep';
end

[num_voxel, num_dim] = size(voxel_sub);
assert(num_dim == 3);
assert(num_voxel >= 2, 'More than 1 pixels are required for computing the orientation');

if num_voxel == 2 && strcmp(method, 'svd')
    warning('Computing orientation using SVD is ill-defined for 2 voxels. Use end to end vector instead.');
    method = 'ep2ep';
end

switch method
    case 'ep2ep'
        ori_vec = voxel_sub(end, :) - voxel_sub(1, :);
        assert(any(ori_vec ~= 0), 'Zero orientation vector');
        ori_vec = ori_vec ./ sqrt(ori_vec * ori_vec');
    case 'svd'        
        voxel_sub= voxel_sub - mean(voxel_sub, 1);
        cov_mat = (voxel_sub' * voxel_sub) ./ (num_voxel - 1);
        [svd_vec, svd_val, ~] = svd(cov_mat);
        ori_vec = svd_vec(:, 1);
        if nargout > 1
            varargout{1} = svd_val(1) ./ trace(svd_val);
        end
    otherwise
        error('Unrecognized method');        
end
end