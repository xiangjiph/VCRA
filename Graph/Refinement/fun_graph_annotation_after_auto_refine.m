function [vessel_graph, record] = fun_graph_annotation_after_auto_refine(vessel_graph, vessel_image, vessel_mask_dt, opt)
%% Initialize parameters
if vessel_graph.link.num_cc == 0
    return;
end
image_size = vessel_graph.num.mask_size;
max_rc_int = opt.mask2graph.max_rc_int;
%% Default value
if isfield(opt, 'internal_offset') && isscalar(opt.internal_offset)
    internal_offset = opt.internal_offset;
else
    internal_offset = 16;
end
if isfield(opt, 'annotation_on_Q') && isscalar(opt.annotation_on_Q)
    annotation_onQ = opt.annotation_on_Q;
else
    annotation_onQ = false;
end

if ~isfield(opt, 'skel_recentering_Q')
    opt.skel_recentering_Q = true;
end
% Record refinement process
if ~isfield(opt, 'record') || ~isfield(opt.record, 'num_self_loop')
    record = struct;
    [record.num_self_loop, record.num_bilink_loop, record.num_short_link_ep1, ...
        record.num_link_ep1, record.num_dim_link_create_ep_if_deleted, record.num_dim_link_change_length, ...
        record.num_dim_link] = deal(0);
    
    [record.dim_short_link_kept_ind, record.invalid_linker_ep_pair, ...
        record.added_linker_ind_o_m, record.added_linker_ind_w_ep, record.int_link, ...
        record.dim_short_link_not_candidate_ind] = deal([]);    
    
    record.finish_refinement_Q = false;
else
    record = opt.record;
end
% Classifier
if ~opt.annotation_on_Q
    c_link_ep1 = opt.classifier.link_ep1;
    c_dim_short_link = opt.classifier.dim_short_link;
    c_linker = opt.classifier.linker;
elseif isfield(opt, 'record') && isfield(opt.record, 'TD')
    TD = opt.record.TD;
else
    TD = [];
end
%% Remove internal link with 1 endpoint
int_link_0 = fun_graph_get_free_link(vessel_graph, internal_offset, true);
if int_link_0.ep1.num_cc > 0
    tmp_ep1_features = fun_graph_get_link_w_1_ep_features(int_link_0.ep1.link_cc_ind, vessel_image, vessel_mask_dt);
    % Classification
    if ~annotation_onQ
        tmp_feature = table2array(tmp_ep1_features(:, c_link_ep1.used_feature_name));
        tmp_ep1_link_removed_Q = c_link_ep1.classifier.predict(tmp_feature);
    else
        % For visualization
        vessel_skl_rc = false(image_size);
        vessel_skl_rc(cat(1, vessel_graph.link.pos_ind, vessel_graph.node.pos_ind, ...
            vessel_graph.isopoint.pos_ind)) = true;
        vessel_mask = vessel_mask_dt > 0;
        
        tmp_check_target = 'link_ep1';
        tmp_annotation_name = sprintf('App_check_%s_result', tmp_check_target);
        tmp = App_label_link(vessel_image, vessel_mask, vessel_skl_rc, int_link_0.ep1.link_cc_ind, tmp_annotation_name);
        waitfor(tmp, 'finishQ', true); 
        try
            annotate_link_ep1_label_1 = tmp.result;
        catch
            annotate_link_ep1_label_1 = eval(tmp_annotation_name);
        end
        disp('Annotation result recorded')
        clear(tmp_annotation_name)
        tmp_ep1_link_removed_Q = (annotate_link_ep1_label_1.to_removeQ | ...
            annotate_link_ep1_label_1.not_sureQ | annotate_link_ep1_label_1.artefactQ) & ...
            (~annotate_link_ep1_label_1.not_sure_but_keepQ & ...
            ~annotate_link_ep1_label_1.artefact_but_keepQ);
        TD = fun_learning_save_annotation_result(tmp_check_target, annotate_link_ep1_label_1, tmp_ep1_features, int_link_0.ep1.link_cc_ind, TD, false, true);    
    end
    tmp_ep1_link_removed_label = int_link_0.ep1.link_label(tmp_ep1_link_removed_Q);
    tmp_kept_link_cc_ind = int_link_0.ep1.link_cc_ind(~tmp_ep1_link_removed_Q);
    % Remove links with endpoints
    vessel_graph = fun_graph_pruning_by_link_label(vessel_graph, tmp_ep1_link_removed_label);
    % The new labels of links with endpoint that need to be kept
    record.num_link_ep1 = record.num_link_ep1 + numel(tmp_ep1_link_removed_label);
    record.valid_link_ep1_ind = tmp_kept_link_cc_ind;
else
%     kept_link_ep1_label = [];    
end
%% Delete erroneous link due to closeness of capillaries and low z resolution 
% Find links that are: 
% 1. Merge error - due to low z-resolution
%   - Normally short, and the sum of distnace transform at two ends are 
% Candidate pre-selection 
link_near_boundary_str = fun_graph_get_boundary_info_str(vessel_graph, internal_offset);
compute_feature_Q = link_near_boundary_str.link_not_near_boundary_Q & vessel_graph.link.num_voxel_per_cc <= opt.select_dim_short_link.max_length;
if isfield(record, 'dim_short_link_not_candidate_ind') && ~isempty(record.dim_short_link_not_candidate_ind)
    tmp_label_str = fun_graph_get_label_of_cc(vessel_graph.link.map_ind_2_label, record.dim_short_link_not_candidate_ind, false);
    compute_feature_Q(tmp_label_str.unique_cc_label_list) = false;
end
% If exist links that are previously labeled as valid, exclude them here.
if ~isempty(record.dim_short_link_kept_ind)
    valid_dim_short_link_label_str = fun_graph_get_label_of_cc(vessel_graph.link.map_ind_2_label, record.dim_short_link_kept_ind, false);
    compute_feature_Q(valid_dim_short_link_label_str.cc_label_list) = false;
%     dim_short_link_label = setdiff(dim_short_link_label, valid_dim_short_link_label_str.cc_label_list);
end
if ~isempty(record.added_linker_ind_o_m)
    added_linker_link_label_str = fun_graph_get_label_of_cc(vessel_graph.link.map_ind_2_label, record.added_linker_ind_o_m, false);
    compute_feature_Q(added_linker_link_label_str.cc_label_list) = false;
%     dim_short_link_label = setdiff(dim_short_link_label, added_linker_link_label_str.cc_label_list);
end
link_features = fun_graph_get_link_features(vessel_graph.link.cc_ind(compute_feature_Q), vessel_image, vessel_mask_dt);
% Select links candidate for checking
is_self_loop_Q = (vessel_graph.link.connected_node_label(compute_feature_Q,1) == vessel_graph.link.connected_node_label(compute_feature_Q,2));
no_endpoint_Q = all(vessel_graph.link.connected_node_label(compute_feature_Q, :), 2);
longer_than_1_Q = link_features.length > 1;

dim_short_link_Q =  link_features.int_min < opt.select_dim_short_link.min_int_max & ... 
    link_features.dt_ep_sum_2_ep_dist > opt.select_dim_short_link.dt_ep_sum_2_ep_dist_min & ...
    link_features.length <= opt.select_dim_short_link.max_length; 

dim_short_link_Q = (~is_self_loop_Q) & no_endpoint_Q & longer_than_1_Q & ...
    dim_short_link_Q;
% Record link whose features do not need to be computed in next round of refinement
compute_feature_label = find(compute_feature_Q);
not_compute_feature_Q = (~compute_feature_Q);
not_compute_feature_Q(compute_feature_label(~dim_short_link_Q)) = true;
record.dim_short_link_not_candidate_ind = vessel_graph.link.cc_ind(not_compute_feature_Q);
%
dim_short_link_label = compute_feature_label(dim_short_link_Q);

dim_short_link_length = link_features.length(dim_short_link_Q);

[~, tmp_idx] = sort(dim_short_link_length, 'descend');
dim_short_link_label = dim_short_link_label(tmp_idx);

dim_short_link_sub_list_ind = find(dim_short_link_Q);
dim_short_link_sub_list_ind = dim_short_link_sub_list_ind(tmp_idx);
% Get loop feature for pre-selected links 
tmp_vessel_graph = fun_analysis_get_connectivity_graph(vessel_graph);
tmp_dim_short_link_loop = fun_analysis_get_loops_in_graph_by_link_label(tmp_vessel_graph, dim_short_link_label, 'euclidean');

dim_short_link_features = link_features(dim_short_link_sub_list_ind, :);
assert(numel(dim_short_link_label) == numel(tmp_dim_short_link_loop.link_label), 'Number of input links are not the same as the output links. Debug');
dim_short_link_features.shortest_loop_length = tmp_dim_short_link_loop.loop_length;
dim_short_link_features.shortest_loop_geodesic_length = tmp_dim_short_link_loop.loop_geodesic_length;
dim_short_link_features.length_ratio_in_shortest_loop = dim_short_link_features.length ./ dim_short_link_features.shortest_loop_length;
dim_short_link_features.ep2ep_vec_wrt_z = dim_short_link_features.ep1_to_ep2_direction_vec(:,3);
if ~isempty(dim_short_link_label)
    if annotation_onQ
        % Update skeleton for visualization
        vessel_skl_rc = false(image_size);
        vessel_skl_rc(vessel_graph.link.pos_ind) = true;
        vessel_skl_rc(vessel_graph.node.pos_ind) = true;
        vessel_skl_rc(vessel_graph.isopoint.pos_ind) = true;

        tmp_check_cc_ind = vessel_graph.link.cc_ind(dim_short_link_label);
        tmp_check_target = 'dim_short_link';
        tmp_annotation_name = sprintf('App_check_%s_result', tmp_check_target);
        tmp = App_label_link(vessel_image, vessel_mask, vessel_skl_rc, tmp_check_cc_ind, tmp_annotation_name);
        waitfor(tmp, 'finishQ', true);
        try
            annotate_link_dim_short_label_1 = tmp.result;
        catch
            annotate_link_dim_short_label_1 = eval(tmp_annotation_name);
        end
        disp('Annotation result recorded')
        clear(tmp_annotation_name)
        dim_short_link_to_remove_Q = annotate_link_dim_short_label_1.to_removeQ | annotate_link_dim_short_label_1.artefactQ | ...
            annotate_link_dim_short_label_1.not_sureQ & ~annotate_link_dim_short_label_1.not_sure_but_keepQ & ...
            ~annotate_link_dim_short_label_1.artefact_but_keepQ;
        TD = fun_learning_save_annotation_result(tmp_check_target, annotate_link_dim_short_label_1, dim_short_link_features, tmp_check_cc_ind, TD, false, true);
    else
        tmp_feature = table2array(dim_short_link_features(:, c_dim_short_link.used_feature_name));
        dim_short_link_to_remove_Q = c_dim_short_link.classifier.predict(tmp_feature);
    end
    
    dim_short_link_label_to_remove = dim_short_link_label(dim_short_link_to_remove_Q);
    % Record dim_short_link_label_to_keep for next iteration of refinement
    record.dim_short_link_kept_ind = cat(1, record.dim_short_link_kept_ind, vessel_graph.link.cc_ind(dim_short_link_label(~dim_short_link_to_remove_Q)));
    % Remove links without creating endpoint 
    [vessel_graph, tmp_undeleted_cc] = fun_graph_delete_dim_short_link(vessel_graph, ...
            dim_short_link_label_to_remove);    
    
    record.num_dim_link_create_ep_if_deleted = record.num_dim_link_create_ep_if_deleted + nnz(tmp_undeleted_cc.will_create_endpointQ);
    record.num_dim_link_change_length = record.num_dim_link_change_length + nnz(tmp_undeleted_cc.change_lengthQ);
    record.num_dim_link = record.num_dim_link + numel(dim_short_link_label_to_remove) - record.num_dim_link_change_length - ...
        record.num_dim_link_create_ep_if_deleted;
end
% Remaining internal link 
int_link_1 = fun_graph_get_free_link(vessel_graph, internal_offset, true);
if int_link_1.num.ep > 0
%% Find linkers
    if isempty(record.added_linker_ind_o_m)
        tmp_ind = cat(1, vessel_graph.endpoint.pos_ind, vessel_graph.link.pos_ind, ...
            vessel_graph.node.pos_ind, vessel_graph.isopoint.pos_ind);
    else
        % Delete isolated voxels if not in the first round of iteration
        tmp_ind = cat(1, vessel_graph.endpoint.pos_ind, vessel_graph.link.pos_ind, ...
            vessel_graph.node.pos_ind);
    end
    vessel_skl = zeros(image_size, 'single');
    vessel_skl(tmp_ind) = max(1, vessel_mask_dt(tmp_ind));
    vessel_mask_1 = vessel_mask_dt | vessel_skl;
    
    [linker_str, th_rlx_ep_pair, ~] = fun_graph_get_linker_for_int_link(vessel_image, vessel_mask_1, vessel_skl, int_link_1);
    nearest_ep_pair = fun_graph_find_nearest_endpoint_pair(int_link_1);
    nearest_ep_pair_ind = cat(1, nearest_ep_pair.ep1_ep1.pair_ind, ...
        nearest_ep_pair.ep1_ep2.pair_ind, nearest_ep_pair.ep2_ep2.pair_ind);
    [ep2ep_pair_ind, ~] = fun_setdiff_row_unordered(nearest_ep_pair_ind, th_rlx_ep_pair);
    [ep2ep_linker_str, ~] = fun_graph_get_linker_for_endpoint_pairs(vessel_image, vessel_mask_1, vessel_skl, ep2ep_pair_ind);
    linker_str = cat(2, linker_str, ep2ep_linker_str);
    if ~isempty(linker_str)
        linker_features = fun_graph_get_linker_features(linker_str, vessel_graph);
        %% Annotate linker    
        if annotation_onQ
            tmp_check_cc_ind = {linker_str.link_ind_w_ep};
            tmp_check_target = 'linker';
            tmp_annotation_name = sprintf('App_check_%s_result', tmp_check_target);
            tmp = App_label_link(vessel_image, vessel_mask_dt > 0, vessel_skl > 0, tmp_check_cc_ind, tmp_annotation_name);
            waitfor(tmp, 'finishQ', true);
            try
                annotate_linker_label = tmp.result;
            catch
                annotate_linker_label = eval(tmp_annotation_name);
            end
            disp('Annotation result recorded')
            clear(tmp_annotation_name)
            valid_linker_Q = (~annotate_linker_label.artefactQ & ~annotate_linker_label.to_removeQ & ...
                ~annotate_linker_label.not_sureQ ) | annotate_linker_label.not_sure_but_keepQ | ...
                annotate_linker_label.artefact_but_keepQ;
            % Save annotation result as training set
            TD = fun_learning_save_annotation_result(tmp_check_target, annotate_linker_label, linker_features, linker_str, TD, false, true);
        else
            valid_linker_Q = ~c_linker.classifier.predict(table2array(linker_features(:, c_linker.used_feature_name)));
        end
        %% Add linkers
        valid_linker_str = linker_str(valid_linker_Q);
        % Record for next iteration of refinment
        invalid_ep_ind_1 = [linker_str(~valid_linker_Q).ep_1_ind]';
        invalid_ep_ind_2 = [linker_str(~valid_linker_Q).ep_2_ind]';
        invalid_linker_ep_pair = cat(2, invalid_ep_ind_1, invalid_ep_ind_2);

        % Endpoint pairs that are classified to be unable to find valid linkers
        record.invalid_linker_ep_pair = cat(1, record.invalid_linker_ep_pair, invalid_linker_ep_pair);
        % Added linker ind
        record.added_linker_ind_o_m = cat(1, record.added_linker_ind_o_m, {valid_linker_str.link_ind_o_m}');
        record.added_linker_ind_w_ep = cat(1, record.added_linker_ind_w_ep, {valid_linker_str.link_ind_w_ep}');

        % Find the linker that have close by endpoints or common endpoints
        valid_linker_str = fun_graph_nearby_linker_selection(valid_linker_str, 10);

        % Add linker as voxel list to the skeleton.
        vessel_skl = vessel_skl > 0;
        if ~isempty(valid_linker_str)
            vessel_skl(cat(1, valid_linker_str.link_ind_w_ep)) = true;
            vessel_graph = fun_skeleton_to_graph(vessel_skl);
        end
    end
end
if opt.annotation_on_Q
    record.TD = TD;
end
%% Deal with the remaining internal endpoints
% int_link_2 = fun_graph_get_free_link(vessel_graph, internal_offset);
% If the link with 2 endpoints all connect to the wrong linkers, and thus
% they are not changed by the linker, remove it
% Record for next iteration of refinement
int_link_2 = fun_graph_get_free_link(vessel_graph, internal_offset);
% Check if the current int_link_2 is the same as the previous one recorded.
if ~isempty(record.int_link)
   ep_ind_intersect = intersect(int_link_2.ep_ind, record.int_link.ep_ind);
   if numel(ep_ind_intersect) == numel(record.int_link.ep_ind)
       record.finish_refinement_Q = true;
   else
       record.int_link = int_link_2;
   end
else
    record.int_link = int_link_2;
end
%% Debug
% assignin('base', 'vessel_graph_1', vessel_graph);
% assignin('base', 'round_1_record', record);
end