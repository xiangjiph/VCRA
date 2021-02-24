function region_graph = fun_graph_construct_graph_from_skl_list(skl_cell)
% fun_graph_construct_graph_in_region_by_grid_sub construct the vascular
% graph in the coordinate of the entire dataset at 1 um voxel resolution. 
% 
is_valid_Q = ~cellfun(@isempty, skl_cell);
skl_cell = skl_cell(is_valid_Q);
num_bbox = nnz(is_valid_Q);

skel_global_ind = cell(num_bbox, 1);
skel_r = cell(num_bbox, 1);
data_size = [];
for iter_bbox = 1 : num_bbox
   try 
       tmp_skl_str = skl_cell{iter_bbox};
       
       if isempty(data_size) && isfield(tmp_skl_str, 'dataset_size')
           data_size = tmp_skl_str.dataset_size;
       end
       tmp_ind = tmp_skl_str.ind;
       
       if ~isempty(tmp_ind)
           tmp_sub = fun_ind2sub(tmp_skl_str.block_size, tmp_ind);
           tmp_sub_g = bsxfun(@plus, tmp_sub, tmp_skl_str.global_bbox_mmll(1:3) - 1);
           tmp_ind_g = sub2ind(tmp_skl_str.dataset_size, tmp_sub_g(:, 1), ...
               tmp_sub_g(:, 2), tmp_sub_g(:, 3));       
           skel_global_ind{iter_bbox} = tmp_ind_g;
           skel_r{iter_bbox} = tmp_skl_str.r;
       end
   catch ME
       fprintf('Failed to load the skeleton of block (%d, %d, %d)\n', tmp_idx_1, tmp_idx_2, tmp_idx_3);
       fun_rethrow_error_message_without_error(ME);
   end      
end
tmp_tic = tic;
skel_global_ind = cat(1, skel_global_ind{:});
skel_r = cat(1, skel_r{:});
[skel_global_ind, unique_idx, ~] = unique(skel_global_ind, 'stable');
skel_r = skel_r(unique_idx);
region_graph = fun_skeleton_to_graph(skel_global_ind, data_size);
% An alternative way is to use fun_graph_add_radius
region_graph.radius = sparse(skel_global_ind, 1, double(skel_r), prod(data_size), 1);
fprintf('Finishing reconstructing the graph. Elapsed time is %f seconds\n', toc(tmp_tic));
end