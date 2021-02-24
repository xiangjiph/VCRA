function [ep2ep_linker_str, ep2ep_linker_found_idx] = fun_graph_get_linker_for_endpoint_pairs(vessel_image, vessel_mask, vessel_skl, ep_pair_ind)
% fun_graph_get_linker_for_endpoint_pairs computes the linker for given
% endpoint pairs by threshold relaxation 
%
%
% Parameters
connectivity = 26;

ep2ep_linker_str = [];
ep2ep_linker_found_idx = [];
image_size = size(vessel_image);
if ~isempty(ep_pair_ind)
    if iscolumn(ep_pair_ind)
        ep_pair_ind = ep_pair_ind';
    end
    sub_1_list = fun_ind2sub(image_size, ep_pair_ind(:,1))';
    sub_2_list = fun_ind2sub(image_size, ep_pair_ind(:,2))';
else
    return;
end
clearvars ep2ep_linker_str
num_pair = size(sub_2_list,2);
if num_pair 
    for iter_pair = num_pair:-1:1
        ep2ep_linker_str(iter_pair) = fun_graph_get_linker_ep2ep_th_rlx(vessel_image, vessel_mask, vessel_skl, ...
            sub_1_list(:, iter_pair), sub_2_list(:, iter_pair), connectivity);
    end
    ep2ep_linker_found_idx = find([ep2ep_linker_str.foundQ]);
    ep2ep_linker_str = ep2ep_linker_str(ep2ep_linker_found_idx);
end
end