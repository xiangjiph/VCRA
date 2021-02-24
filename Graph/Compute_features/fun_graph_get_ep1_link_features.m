function lf = fun_graph_get_ep1_link_features(link_label_list, vessel_graph, vessel_image, vessel_mask_dt)

% Parameters
local_bbox_expansion = 10;
% Initialization
lf = struct;
lf.num_l = numel(link_label_list);
lf.num_voxel = vessel_graph.link.num_voxel_per_cc(link_label_list);
lf.length = zeros(lf.num_l, 1);
lf.ep_dist = ones(lf.num_l, 1);
lf.dt_max = zeros(lf.num_l, 1);
lf.dt_min = zeros(lf.num_l, 1);
lf.dt_mean = zeros(lf.num_l, 1);
lf.dt_std = zeros(lf.num_l, 1);
lf.int_mean = zeros(lf.num_l, 1);
lf.int_median = zeros(lf.num_l, 1);
lf.int_std = zeros(lf.num_l, 1);
lf.node_int = zeros(lf.num_l, 1);
% Information about the links that share the same node as this link
lf.nl_int_median = zeros(lf.num_l, 1);
lf.nl_int_std = zeros(lf.num_l, 1);
lf.nl_max_length = zeros(lf.num_l, 1);
lf.nl_dt_mean = zeros(lf.num_l, 1);
lf.nl_dt_std = zeros(lf.num_l, 1);
% Information about the background statistics
lf.bg_mean = zeros(lf.num_l, 1);
lf.bg_std = zeros(lf.num_l, 1);
%% 
for iter_link = 1 : lf.num_l
    tmp_label = link_label_list(iter_link);
    tmp_ind = vessel_graph.link.cc_ind{tmp_label};
    tmp_sub = fun_ind2sub(vessel_graph.num.mask_size, tmp_ind);
    lf.length(iter_link) = fun_graph_ind_to_length(tmp_sub, 1);
    if size(tmp_sub,1) > 1
        lf.ep_dist(iter_link) = sqrt(sum((tmp_sub(1, :) - tmp_sub(end,:)).^2));
    end
    tmp_dt = vessel_mask_dt(tmp_ind);
    tmp_int = vessel_image(tmp_ind);
    lf.dt_max(iter_link) = max(tmp_dt);
    lf.dt_min(iter_link) = min(tmp_dt);
    lf.dt_mean(iter_link) = mean(tmp_dt);
    lf.dt_std(iter_link) = std(tmp_dt);
    lf.int_mean(iter_link) = mean(tmp_int);
    lf.int_median(iter_link) = median(tmp_int);
    lf.int_std(iter_link) = std(single(tmp_int));
    % Neighboring link information
    tmp_node_label = vessel_graph.link.connected_node_label(tmp_label,:);
        % For debug
    assert(nnz(tmp_node_label) == 1, 'The link connected to two nodes!');
    tmp_node_label = tmp_node_label(1);
    tmp_neighbor_link_label = vessel_graph.node.connected_link_label{tmp_node_label};
    tmp_neighbor_link_label = tmp_neighbor_link_label(tmp_neighbor_link_label~=tmp_label);
    assert(length(tmp_neighbor_link_label) >= 2, 'The number of neighboring link is less than 2');
    % Neighboring link 
    tmp_nl_cc = vessel_graph.link.cc_ind(tmp_neighbor_link_label);
    tmp_nl_ind = cat(1, tmp_nl_cc{:});
    lf.nl_max_length(iter_link) = max(cellfun(@length, tmp_nl_cc));
    tmp_nl_int = vessel_image(tmp_nl_ind);
    tmp_nl_dt = vessel_mask_dt(tmp_nl_ind);
    lf.nl_int_median(iter_link) = median(tmp_nl_int);
    lf.nl_int_std(iter_link) = std(single(tmp_nl_int));
    lf.nl_dt_mean(iter_link) = mean(tmp_nl_dt);
    lf.nl_dt_std(iter_link) = std(tmp_nl_dt);
    % Background information
    tmp_min = max([1,1,1],  min(tmp_sub) - local_bbox_expansion);
    tmp_max = min(vessel_graph.num.mask_size, max(tmp_sub) + local_bbox_expansion);
    tmp_mmll = [tmp_min, tmp_max - tmp_min + 1];
    tmp_mask = crop_bbox3(vessel_mask_dt, tmp_mmll, 'default') > 0;
    tmp_image = crop_bbox3(vessel_image, tmp_mmll, 'default');
    tmp_bg_int = tmp_image(~tmp_mask);
    lf.bg_mean(iter_link) = mean(tmp_bg_int);
    lf.bg_std(iter_link) = std(single(tmp_bg_int));
end
lf.straightness = lf.ep_dist ./ lf.length;
lf.int_SNR = (lf.int_median - lf.bg_mean) ./ lf.int_std;
lf.dt_std_n = lf.dt_std ./ lf.dt_mean;
end