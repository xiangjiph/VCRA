function ep_pair = fun_graph_get_nearest_point_pair_1(ep_sub, search_range_th)




if nargin < 2
    search_range_th = inf;
end

dist_p2p = squareform(pdist(ep_sub));
[mutual_nn_1_idx, mutual_nn_2_idx, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_p2p, true);
tmp_Q = mutual_nn_dist <= search_range_th;
mutual_nn_1_idx = mutual_nn_1_idx(tmp_Q);
mutual_nn_2_idx = mutual_nn_2_idx(tmp_Q);
mutual_nn_dist = mutual_nn_dist(tmp_Q);

[~, tmp_unique_idx ] = unique(sort(cat(2, mutual_nn_1_idx, mutual_nn_2_idx), 2, 'ascend'), 'row', 'stable');

ep_pair.pair_idx = cat(2, mutual_nn_1_idx(tmp_unique_idx), mutual_nn_2_idx(tmp_unique_idx));
ep_pair.pair_sub_1 = ep_sub(ep_pair.pair_idx(:,1), :);
ep_pair.pair_sub_2 = ep_sub(ep_pair.pair_idx(:,2), :);
ep_pair.pair_dist = mutual_nn_dist(tmp_unique_idx);
end