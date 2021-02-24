function pair_str = fun_graph_get_nearby_point_pairs(point_sub, search_range)

% search_range = 50;% um

% Compute pair-wise distance
dist_ep1_to_ep1 = squareform(pdist(point_sub));
% Find nearby points
nearby_Q = (dist_ep1_to_ep1 <= search_range) & (dist_ep1_to_ep1 ~= 0); % Remove self-distance(0);
% Index of point with nearby points 
p2p_idx = find(any(nearby_Q, 2));
if isempty(p2p_idx)
    pair_str.idx = [];
    return;
end
p2p_paired_idx_cell = cell(numel(p2p_idx),1);
p2p_paired_dist_cell = cell(numel(p2p_idx),1);
for iter_ep = 1 : numel(p2p_idx)
    p2p_paired_idx_cell{iter_ep} = find(nearby_Q(p2p_idx(iter_ep), :))';
    p2p_paired_dist_cell{iter_ep} = dist_ep1_to_ep1(p2p_idx(iter_ep),p2p_paired_idx_cell{iter_ep})';
end
p2p_num_pair = cellfun(@numel, p2p_paired_idx_cell);
p2p_idx_pair = cat(2, repelem(p2p_idx, p2p_num_pair), cat(1, p2p_paired_idx_cell{:}));

p2p_dist = cat(1, p2p_paired_dist_cell{:});
[p2p_idx_pair_unique, unique_idx]= unique(sort(p2p_idx_pair, 2, 'ascend'), 'row', 'stable');
pair_str.idx = p2p_idx_pair_unique;
pair_str.dist = p2p_dist(unique_idx);
pair_str.sub_1 = point_sub(p2p_idx_pair_unique(:,1), :);
pair_str.sub_2 = point_sub(p2p_idx_pair_unique(:,2), :);
end