function varargout = fun_graph_point_pairs_is_significantly_unique(sub_1_list, sub_2_list, dist_th)
% When finding the linker for gaps, we find all the possible linker first,
% then do the selection. However, sometimes the pair of points can be very
% close to another pair of points. This function finds the nearby point
% pairs and remove one appears later in the list. 
% Input: 
%   sub_1_list: N-by-3 numerical array, subscript of the points
%   sub_2_list: N-by-3 numerical array, subscript of the points
%   dist_th: Euclidean distance threshold
% Output: 
%   keptQ: N-by-1 logicla array, true if the point pair has not nearby
%   neighbor pairs
% dist_th = 20 is too large even for selecting linkers
if nargin < 3
    dist_th = 3;
end

if size(sub_1_list, 1) == 3
    warning('Expected input point list to be N-by-3. Transpose automatically');
    sub_1_list = sub_1_list';
end
if size(sub_2_list, 1) == 3
    warning('Expected input point list to be N-by-3. Transpose automatically');
    sub_2_list = sub_2_list';
end
num_pair = size(sub_1_list, 1);
% Nearby endpoints
tmp_sub_list = cat(1, sub_1_list, sub_2_list);
dist_ep_ep = squareform(pdist(tmp_sub_list));
for iter_1 = 1 : num_pair
    % Remove diagonal term
    dist_ep_ep(iter_1, iter_1) = inf;
    dist_ep_ep(iter_1 + num_pair, iter_1 + num_pair) = inf;
    % Remove paired points
    dist_ep_ep(iter_1, iter_1 + num_pair) = inf;
    dist_ep_ep(iter_1 + num_pair, iter_1) = inf;
end
% If both points in the paired voxel are very closed to each other, merge
% them. 
[tmp_idx1, tmp_idx2] = find(dist_ep_ep < dist_th);
tmp_idx1_idx2 = cat(2, tmp_idx1, tmp_idx2);
tmp_idx1_idx2 = unique(sort(tmp_idx1_idx2, 2, 'ascend'), 'rows', 'stable');
tmp_num_pair_to_test = nnz(tmp_idx1_idx2(:,1) <= num_pair);
kept_pair_Q = true(num_pair, 1);
for iter_1 = 1 : tmp_num_pair_to_test
    tmp_idx_1_1 = tmp_idx1_idx2(iter_1,1);
    if kept_pair_Q(tmp_idx_1_1)
        tmp_idx_2_1 = tmp_idx1_idx2(iter_1,2);
        if tmp_idx_2_1 > num_pair
            tmp_idx_2_2 = tmp_idx_2_1;
            tmp_idx_2_1 = tmp_idx_2_2 - num_pair;
        else
            tmp_idx_2_2 = tmp_idx_2_1 + num_pair;
        end
        tmp_idx_1_2 = tmp_idx_1_1 + num_pair;    
        if dist_ep_ep(tmp_idx_1_2, tmp_idx_2_2) < dist_th
            fprintf('Pair (%d, %d) and (%d, %d) are duplicated\n', tmp_idx_1_1, tmp_idx_1_2, tmp_idx_2_1, tmp_idx_2_2);
            kept_pair_Q(tmp_idx_2_1) = false;
        end
    end
end

if nargout == 1
    varargout{1} = kept_pair_Q;
elseif nargout == 2
    varargout{1} = sub_1_list(kept_pair_Q, :);
    varargout{2} = sub_2_list(kept_pair_Q, :);
elseif nargout == 3
    varargout{1} = sub_1_list(kept_pair_Q, :);
    varargout{2} = sub_2_list(kept_pair_Q, :);
    varargout{3} = kept_pair_Q;
end

end