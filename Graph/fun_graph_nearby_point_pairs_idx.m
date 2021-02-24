function kept_pair_idx_unique = fun_graph_nearby_point_pairs_idx(sub_1_list, sub_2_list, dist_th)
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

if all(size(sub_1_list) == [3, 1])
    warning('Expected input point list to be N-by-3. Transpose automatically');
    sub_1_list = sub_1_list';
end
if all(size(sub_1_list) == [3, 1])
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
kept_pair_idx = zeros(2, tmp_num_pair_to_test);
for iter_1 = 1 : tmp_num_pair_to_test
    tmp_idx_1_1 = tmp_idx1_idx2(iter_1,1);    
        tmp_idx_2_1 = tmp_idx1_idx2(iter_1,2);
        if tmp_idx_2_1 > num_pair
            tmp_idx_2_2 = tmp_idx_2_1;
            tmp_idx_2_1 = tmp_idx_2_2 - num_pair;
        else
            tmp_idx_2_2 = tmp_idx_2_1 + num_pair;
        end
        tmp_idx_1_2 = tmp_idx_1_1 + num_pair;    
        if dist_ep_ep(tmp_idx_1_2, tmp_idx_2_2) < dist_th
%             fprintf('Pair (%d, %d) and (%d, %d) are duplicated\n', tmp_idx_1_1, tmp_idx_1_2, tmp_idx_2_1, tmp_idx_2_2);
            kept_pair_idx(:, iter_1) = [tmp_idx_1_1;tmp_idx_2_1];
        end
end

kept_pair_idx_unique = kept_pair_idx(:, all(kept_pair_idx>0, 1));
kept_pair_idx_unique = unique(sort(kept_pair_idx_unique, 1, 'ascend')', 'rows', 'stable');
end