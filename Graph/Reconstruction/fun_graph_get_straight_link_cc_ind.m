function link_cc_straight = fun_graph_get_straight_link_cc_ind(link_cc, mask_size)
% fun_graph_get_straight_link_cc_ind compute the approxiamte straight line
% voxel indices for real skeleton link connected components. 
% Input: 
%   link_cc: cell array, each cell contains the linear indices of the
%   location of the link voxels in the mask 
%   mask_size: 3-by-1 numerical vector, the size of the mask 
% Output: 
%   link_cc_straight: cell array, each cell contains the linear indices of
%   the locaiton of the approximate straight line in the mask 
% 
% Implemented by Xiang Ji on 02/23/2019
% 
% Replace the link with a straight line joining the two node. Do the
% distance transform and compare the tissue to vessel statistics

num_cc = numel(link_cc);
link_cc_straight = cell(num_cc, 1);
for iter_cc = 1 : num_cc
    tmp_ind = link_cc{iter_cc};
    if ~isempty(tmp_ind)
        if numel(tmp_ind) <= 2
            link_cc_straight{iter_cc} = tmp_ind;
        else
            tmp_sub = fun_ind2sub(mask_size, tmp_ind);
            tmp_ep1_sub = tmp_sub(1, :);
            tmp_ep2_sub = tmp_sub(end, :);
            % The following approach is just the approximation of the line
            % that connect the two ends of the link. The resulting indice
            % list is not a link cc because it probably contains removable
            % points. However, this approach is much faster than the
            % approach that is commented out below.  
            tmp_ep2ep_vec = tmp_ep2_sub - tmp_ep1_sub;
            tmp_dist_ep2ep = sqrt(sum(tmp_ep2ep_vec.^2));
            tmp_ep2ep_vec = tmp_ep2ep_vec ./ tmp_dist_ep2ep;
            tmp_straight_line_sub = round(tmp_ep1_sub +  bsxfun(@times, (0:1:tmp_dist_ep2ep)', tmp_ep2ep_vec));
%             tmp_straight_line_sub = unique(tmp_straight_line_sub, 'rows', 'stable');
            tmp_valid_sub_Q = all(tmp_straight_line_sub > 0, 2) & all(tmp_straight_line_sub < mask_size, 2);
            tmp_straight_line_sub = tmp_straight_line_sub(tmp_valid_sub_Q, :);
            link_cc_straight{iter_cc} = sub2ind(mask_size, tmp_straight_line_sub(:,1), ...
                tmp_straight_line_sub(:,2), tmp_straight_line_sub(:,3));
            
            % Compute in the local coordinate to save time...
%             tmp_sub_min = min(tmp_sub, [], 1);
%             tmp_sub_max = max(tmp_sub, [], 1);
%             tmp_sub_max = max(tmp_sub_max, tmp_sub_min + 2);
%             bbox_size = tmp_sub_max - tmp_sub_min + 1;
%             tmp_ep1_sub_local = tmp_ep1_sub - tmp_sub_min + 1;
%             tmp_ep2_sub_local = tmp_ep2_sub - tmp_sub_min + 1;
%                        
%             tmp_mask = fun_graph_get_p2p_cylinder_mask(bbox_size, tmp_ep1_sub_local, ...
%                 tmp_ep2_sub_local, 1);
%             shortest_path_str = fun_graph_shortest_path_between_points(tmp_mask, tmp_ep1_sub_local, ...
%                 tmp_ep2_sub_local, 26);
%             linker_sub = shortest_path_str.link_sub_w_ep + tmp_sub_min - 1;
%             link_cc_straight{iter_cc} = sub2ind(mask_size, linker_sub(:,1), linker_sub(:,2), linker_sub(:,3));
        end
    end
end
end