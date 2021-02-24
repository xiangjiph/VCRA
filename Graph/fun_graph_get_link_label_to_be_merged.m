function link_label = fun_graph_get_link_label_to_be_merged(cc_ind, mask_dt, max_merge_length, always_merge_max_length)
% fun_graph_get_link_label_to_be_merged computes the geometric features of
% the input link segments and determines the label of the links to be
% merged. 
% Input: 
%    cc_ind: cell of link voxle linear indices
%   mask_dt: 3D single precision array, distance transform of the vessel
%   mask 
%   max_merge_length: link of length greater than this value won't be
%   merged no matter what. 
%   always_merge_max_legth: link of length no greater than this value will be
%   merged no matter what. 
% Ouptut: 
%   link_label: numericla vector, indices of the link to be removed in the
%   input cc_ind
if nargin < 3
    max_merge_length = inf;
    always_merge_max_length = 0;
elseif nargin < 4
    always_merge_max_length = 0;
end
dt_is_sparseQ = issparse(mask_dt);
image_size = size(mask_dt);
lf = struct;
num_l = numel(cc_ind);
if num_l == 0
    link_label = [];
    return;
end
lf.num_voxel = cellfun(@length, cc_ind);
[lf.length, lf.dt_max ] = deal(nan(num_l, 1));
%% Compute single link features
for iter_link = 1 : num_l
    tmp_ind = cc_ind{iter_link};
    tmp_sub = fun_ind2sub(image_size, tmp_ind);
    lf.length(iter_link) = fun_graph_sub_to_length(tmp_sub, 1);
    tmp_num_voxel = lf.num_voxel(iter_link);
    if tmp_num_voxel > 1
        lf.ep2ep_dist(iter_link) = sqrt(sum((tmp_sub(1, :) - tmp_sub(end,:)).^2));
    end
    tmp_dt = mask_dt(tmp_ind);
    if dt_is_sparseQ
        tmp_dt = full(tmp_dt);
    end
    lf.dt_max(iter_link) = max(tmp_dt);
end
link_label = find((lf.length <= min(lf.dt_max, max_merge_length)) | (lf.num_voxel <= always_merge_max_length));
end