function endpoint_feature = fun_graph_get_endpoint_features(ep_sub, ep_vec)
% fun_graph_get_endpoint_features computes the features derived from the
% position and orientation of endpoints and their nearest neighbor
% Input: 
%   ep_sub: N-by-3 numerical array, where N is the number of endpoint and 3
%   is the 3D subscript of the endpoint in the array
%   ep_vec: 3-by-N numerical array, the normalize direction vector of the
%   endpoint, which is computed as the first normalized principle component
%   of the link that the endpoint attached to, which is pointing along the link and
%   ends at the endpoint. 
% Output: 
%   endpoint_feature: structure with fields
%       nearest_ep_dist: distance to the closest endpoint 
%       nearest_ep_idx: index of the closest endpoint in the ep_sub list
%       ep_ep_disp_vec_n: normalized displacement vector between two
%       endpoints, pointing from the one to its nearest neighbor
%       inner_product_epv_epv: inner product of two direction vectors of
%       the endpoint
%       inner_product_epv1_dispv: inner product of the first endpoint
%       direction vector and the endpoint displacement vector
%       inner_product_epv2_dispv: inner product of the second endpoint
%       direction vector and the endpoint displacement vector (flip
%       direction)

dist_ep2ep = squareform(pdist(ep_sub));
num_ep = size(ep_sub, 1);
dist_ep2ep = dist_ep2ep + diag(ones(num_ep,1) * inf);
% Find the nearest endpoint for each endpoint, not necessary mutually
% matched - otherwise, there will be missing features
[endpoint_feature.nearest_ep_dist, endpoint_feature.nearest_ep_idx] = min(dist_ep2ep, [], 2);

% Inner produce between the two endpoint vector: 
tmp_matched_ep_vec = ep_vec(:, endpoint_feature.nearest_ep_idx);
tmp_matched_ep_sub = ep_sub(endpoint_feature.nearest_ep_idx, :);
tmp_ep_ep_disp_vec = tmp_matched_ep_sub - ep_sub;
tmp_ep_ep_disp_vec = tmp_ep_ep_disp_vec ./ vecnorm(tmp_ep_ep_disp_vec, 2, 2);
endpoint_feature.ep_ep_disp_vec_n = tmp_ep_ep_disp_vec;

endpoint_feature.inner_product_epv_epv = sum(ep_vec .* tmp_matched_ep_vec, 1)';
endpoint_feature.inner_product_epv1_dispv = sum(tmp_ep_ep_disp_vec' .* ep_vec, 1)';
endpoint_feature.inner_product_epv2_dispv = sum(- tmp_ep_ep_disp_vec' .* tmp_matched_ep_vec, 1)';
end