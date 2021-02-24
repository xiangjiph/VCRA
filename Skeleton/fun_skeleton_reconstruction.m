function mask = fun_skeleton_reconstruction(skeleton_voxel, skeleton_radius, block_size)

% Classify skeleton voxel
max_r_edge = ceil(max(skeleton_radius));
radius_bin_val = [1, sqrt(2), 2, sqrt(5), 3 : max_r_edge];
radius_bin_edge = [0, movmean(radius_bin_val, 2, 'Endpoints', 'discard'), max_r_edge];

r_class_ind_cell = fun_bin_data_to_idx_list_by_edges(skeleton_radius, radius_bin_edge, true);
assert(sum(cellfun(@numel, r_class_ind_cell)) == numel(skeleton_voxel), 'Number of skeleton voxel in the cell array does not equal the number of input skeleton voxel')
% Binary reconstruction
mask = false(block_size);
for tmp_idx = 1 : numel(r_class_ind_cell)
    tmp_ind = r_class_ind_cell{tmp_idx};
    if ~isempty(tmp_ind)
        tmp_radius = radius_bin_val(tmp_idx);   
        switch tmp_radius
            case radius_bin_val(1)
                tmp_ind = skeleton_voxel(tmp_ind);
            case radius_bin_val(2)
                tmp_ind = fun_skeleton_reconstruction_ind(skeleton_voxel(tmp_ind), block_size, strel('sphere',1).Neighborhood);
            case radius_bin_val(3)
                tmp_ind = fun_skeleton_reconstruction_ind(skeleton_voxel(tmp_ind), block_size, strel('cube', 3).Neighborhood);
            case radius_bin_val(4)
                tmp_ind = fun_skeleton_reconstruction_ind(skeleton_voxel(tmp_ind), block_size, strel('sphere', 2).Neighborhood);
            otherwise
                tmp_ind = fun_skeleton_reconstruction_ind(skeleton_voxel(tmp_ind), block_size, strel('sphere', double(tmp_radius)).Neighborhood);
        end
        mask(tmp_ind) = true;
    end
end
    
end

