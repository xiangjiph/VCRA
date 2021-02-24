function [X_,Y_,rate_,pixshift,nonuniformity] = searchpair_vessel(des_fixed,des_moving_ori,pixshiftinit,iadj,tile_size,matchparams)
%SEACHPAIR Summary of this function goes here
%
% [OUTPUTARGS] = SEACHPAIR(INPUTARGS) Explain usage here
%
% Inputs:
%   descent: Num_descriptor-by-Num_feature numerical array of the tile
%   descadjori: Num_descriptor-by-Num_feature numerical array of the
%   neighboring tile
%   pixshiftinit: initially estimated pixel shift
%   iadj: index of adjacant tile that have relative position shift w.r.t.
%   the current tile (?)
%   dims: tile size
%   matchparams: 
%   
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/11/03 16:02:56 $	$Revision: 0.1 $
% Copyright: HHMI 2016
if isfield(matchparams, 'max_num_desc')
    sample_opt.total_num_descriptor = matchparams.max_num_desc;
else
    sample_opt.total_num_descriptor = 5000;
end

if ~isfield(matchparams, 'viz')
    vis_Q = false;
else
    vis_Q = matchparams.viz;
end
if isfield(matchparams, 'remove_descriptor_by_pdist2')
    sample_opt.rm_descriptor_by_pdist2_Q = matchparams.remove_descriptor_by_pdist2;
else
    sample_opt.rm_descriptor_by_pdist2_Q = true;
end


[X_,Y_] = deal([]);
pixshiftinit = round(pixshiftinit);
pixshift = pixshiftinit;
% Inconsistancy threshold
th_inconsistancy = 0.2;
if ~isfield(matchparams, 'pixshift_search_shift_z')
    pixshift_search_shift_z = [1; -1] * [15, 30, 50];
else
    pixshift_search_shift_z = [1; -1] * matchparams.pixshift_search_shift_z;
end
num_search_option = numel(pixshift_search_shift_z);
flag_stop = false;
iter = 0;
switch iadj
    case 1
        sample_opt.max_disp_pixel_yxz = [15, 10, 5];
    case 2
        sample_opt.max_disp_pixel_yxz = [10, 15, 5];
    case 3
        sample_opt.max_disp_pixel_yxz = [30, 30, 30]; % I am not sure if these numbers are good 
end

R_consistant = zeros(1,50);
R_matched = zeros(1,50);
nonuniformity = zeros(1,numel(pixshift_search_shift_z));
% clear nonuniformity
search_idx = 1;
rate_ = 0;
search_by_brute_force_Q = true;
while ~flag_stop && iter <= num_search_option% run a search
    [X_,Y_, ] = deal([]);
    iter = iter + 1;
%% Shift the pixels according to the stage displacement 
    des_1_sub = des_fixed(:, 1:3);
    des_1_label = des_fixed(:, 4);
    des_1_sub_min = min(des_1_sub, [], 1);
    des_1_sub_max = max(des_1_sub, [], 1);
    
    des_2_sub_shift = bsxfun(@plus, des_moving_ori(:, 1:3), round(pixshift));
    des_2_label = des_moving_ori(:, 4);
    des_2_sub_shift_min = min(des_2_sub_shift, [], 1);
    des_2_sub_shift_max = max(des_2_sub_shift, [], 1);
    
    overlap_bbox_max = min([des_1_sub_max; des_2_sub_shift_max; tile_size],[], 1);
    overlap_bbox_min = max([des_1_sub_min; des_2_sub_shift_min; [1,1,1]],[], 1);
    
    desc_1_selected_Q = all(bsxfun(@ge, des_1_sub, overlap_bbox_min) & bsxfun(@le, des_1_sub, overlap_bbox_max), 2);
    desc_2_selected_Q = all(bsxfun(@ge, des_2_sub_shift, overlap_bbox_min) & bsxfun(@le, des_2_sub_shift, overlap_bbox_max), 2);
    
    des_1_sub = des_1_sub(desc_1_selected_Q, :);
    des_2_sub_shift = des_2_sub_shift(desc_2_selected_Q, :);    
    if isempty(des_1_sub) || isempty(des_2_sub_shift)
        return;
    end
    des_1_label = des_1_label(desc_1_selected_Q);
    des_2_label = des_2_label(desc_2_selected_Q);    
%% Sample the voxels for matching
    [des_1_sub, des_2_sub_shift, des_1_label, des_2_label] = fun_feature_match_sample_descriptor(des_1_sub, ...
        des_2_sub_shift, des_1_label, des_2_label, sample_opt);
    
    tmp_ind = sub2ind(tile_size, des_1_sub(:,1), des_1_sub(:,2), des_1_sub(:,3));
    X_ind_2_label = sparse(tmp_ind, ones(length(tmp_ind),1), des_1_label, prod(tile_size),1);
    
    tmp_ind = sub2ind(tile_size, des_2_sub_shift(:,1), des_2_sub_shift(:,2), des_2_sub_shift(:,3));
    Y_ind_2_label = sparse(tmp_ind, ones(length(tmp_ind),1), des_2_label, prod(tile_size),1);
    
%% Coherent Point drift
    if size(des_1_sub,1)<3 || size(des_2_sub_shift,1)<3% not enough sample to match
        flag_stop = 1;
    else
        % check uniformity of data
        nbins = [2 2];
        edges = {};
        for ii = 2:-1:1%length(dims)%[1 2 3],
            minx = 0;
            maxx = tile_size(ii);
            binwidth = (maxx - minx) / nbins(ii);
            edges{ii} = minx + binwidth*(0:nbins(ii));
        end
        % Get the bivaritive histogram for the location of the feature
        % points. These point are counted in four quarters of the x-y plane
        [accArr] = hist3([des_1_sub(:,1:2);des_2_sub_shift(:,1:2)],'Edges',edges);
        accArr = accArr(1:2,1:2);
        % Check if the four quarter are balanced by checking if the
        % diagonal quarters are larger than the average while the
        % anti-diagonal quaters are less than the average, or vice versa.
        if ~all(sum( accArr>mean(accArr(:)) ) & sum(accArr>mean(accArr(:)),2)' )
            % non uniform over quad-representation
            nonuniformity(iter) = 1;
        else
            nonuniformity(iter) = 0;
        end
%% Match the descriptor - nonregid        
        [rate, X_, Y_, ~] = descriptorMatchforz(des_1_sub,des_2_sub_shift,pixshift,iadj,matchparams);
%% Affine transformation 
% opt.method='affine';
% opt.beta=6;            % the width of Gaussian kernel (smoothness)
% opt.lambda=16;          % regularization weight
% opt.viz=1;              % show every iteration
% opt.outliers=0.9;       % use 0.7 noise weight
% opt.max_it = 50;        % Maxinum number of iteration 
% opt.fgt=0;              % do not use FGT (default)
% opt.tol = 1e-3;         % Error tolorance ( how is it defined? )
% opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
% opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
% matchparams_af.opt = opt;
% matchparams_af.projectionThr = 5;
% matchparams_af.model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
% matchparams_af.debug = false;
% matchparams_af.viz = true;
% tic
% [rate_af, X_af, Y_af, tY_af] = vessel_descriptorMatchforz(des_1_sub, des_2_sub_shift, pixshift, iadj, matchparams_af);
% toc   
        %%
        % Check if the matched points are from the same connected
        % components or not. A good matches should contains many voxel
        % pairs from single connected components
        if isempty(X_) || isempty(Y_)
            return;
        end
%% Matched voxel pair selection 
        X_ind = sub2ind(tile_size, X_(:,1), X_(:,2), X_(:,3));
        Y_ind = sub2ind(tile_size, Y_(:,1) + pixshift(1), Y_(:,2) + pixshift(2), Y_(:,3) + pixshift(3));
        X_label = full(X_ind_2_label(X_ind));
        Y_label = full(Y_ind_2_label(Y_ind));
        % Connected components in the input point set
        X_cc_idx = fun_bin_data_to_idx_list(X_label);
        X_cc_num = numel(X_cc_idx);
        X_cc_size = cellfun(@numel, X_cc_idx);
        % Define inconsistancy as the number of labels / number of voxels
        X_cc_inconsistancy = zeros(X_cc_num, 1);
        for iter1 = 1 : X_cc_num
            X_cc_inconsistancy(iter1) = numel(unique(Y_label(X_cc_idx{iter1})))/X_cc_size(iter1);
        end
        matched_cc_x = X_cc_inconsistancy <= th_inconsistancy;
        
        
        Y_cc_idx = fun_bin_data_to_idx_list(Y_label);
        Y_cc_num = numel(Y_cc_idx);
        Y_cc_size = cellfun(@numel, Y_cc_idx);
        % Define inconsistancy as the number of labels / number of voxels
        Y_cc_inconsistancy = zeros(Y_cc_num, 1);
        for iter1 = 1 : Y_cc_num
            Y_cc_inconsistancy(iter1) = numel(unique(Y_label(Y_cc_idx{iter1})))/Y_cc_size(iter1);
        end
        matched_cc_y = Y_cc_inconsistancy <= th_inconsistancy;
        match_x_Q = false(numel(X_ind),1);
        match_x_Q(cat(2, X_cc_idx{matched_cc_x})) = true;
        
        match_y_Q = false(numel(X_ind),1);
        match_y_Q(cat(2, Y_cc_idx{matched_cc_y})) = true;
        matched_Q = match_x_Q & match_y_Q;
        
        if any(matched_cc_x)
            X_ = X_(matched_Q,:);
            Y_ = Y_(matched_Q,:);
        else
            X_ = [];
            Y_ = [];
        end
        
        if size(X_,1) < 3
            rate = 0;
            consistent_rate = 0; % Too less points matched
        else
            consistent_rate = size(X_,1) / numel(X_ind);
        end
        R_matched(iter) = rate;
        R_consistant(iter) = consistent_rate;
        % The displacement between the matched voxels should be consistant.
        % Remove the outlier pairs whose displacements are significant
        % different from the median displacement of the matched pairs
        disp_X_Y = X_ - Y_;
        disp_X_Y_med = median(disp_X_Y,1);
        disp_X_Y_dev = bsxfun(@minus, disp_X_Y , disp_X_Y_med);
        disp_X_Y_dev_std = std(single(disp_X_Y_dev), 1);
        disp_X_Y_tol = min(15, max(5,disp_X_Y_dev_std * 3));
        disp_kept = all(bsxfun(@le, abs(disp_X_Y_dev), disp_X_Y_tol), 2);
        X_ = X_(disp_kept, :);
        Y_ = Y_(disp_kept, :); 
        if ~isempty(X_) && size(X_,1) > 100 %&& search_by_brute_force_Q
            search_by_brute_force_Q = false;
        else 
            search_by_brute_force_Q = true;
        end
%% Re-estiamte the initial pixel shift - Brute force       
        if matchparams.scan_pixshift_Q && rate < 0.95
            disp('Matching not good enough. Shift the overlapping region and search for pairs again');
            if search_by_brute_force_Q
                if iter == 1
                    X_0 = X_;
                    Y_0 = Y_;
                    pixshift = pixshiftinit;
                    pixshift(iadj) = pixshift(iadj) + pixshift_search_shift_z(1);
                elseif iter == 2
                    X_1 = X_;
                    Y_1 = Y_;
                    pixshift = pixshiftinit;
                    pixshift(iadj) = pixshift(iadj) + pixshift_search_shift_z(2);
                elseif iter == 3
                    if R_consistant(3) > R_consistant(1) && R_consistant(3) > R_consistant(2)
                        % Update the best matched pair
                        X_0 = X_;
                        Y_0 = Y_;
                        search_idx = 4;
                        pixshift = pixshiftinit;
                        pixshift(iadj) = pixshift(iadj) + pixshift_search_shift_z(search_idx);
                    elseif R_consistant(2) > R_consistant(1) && R_consistant(2) > R_consistant(3)
                        X_0 = X_1;
                        Y_0 = Y_1;
                        search_idx = 3;
                        pixshift = pixshiftinit;
                        pixshift(iadj) = pixshift(iadj) + pixshift_search_shift_z(search_idx);
                    else
                        flag_stop = true;
                        X_ = X_0;
                        Y_ = Y_0;
                        rate = R_matched(1);
                        pixshift = median(X_ - Y_);
                    end
                else
                    if R_consistant(iter) > R_consistant(iter - 1) && (search_idx + 2) <= num_search_option
                        search_idx = search_idx + 2;
                        X_0 = X_;
                        Y_0 = Y_;
                        pixshift = pixshiftinit;
                        pixshift(iadj) = pixshift(iadj) + pixshift_search_shift_z(search_idx);
                    elseif R_consistant(iter) > R_consistant(iter - 1) && (search_idx + 2) > num_search_option
                        rate = R_matched(iter-1);
                        flag_stop = true;
                        pixshift = median(X_ - Y_);
                    else
                        X_ = X_0;
                        Y_ = Y_0;
                        pixshift = median(X_ - Y_);
                        flag_stop = true;
                        rate = R_matched(iter-1);
                    end
                end
            else
                if iter == 1
                    X_0 = X_;
                    Y_0 = Y_;
                    pixshift = disp_X_Y_med;
                else
                    if iter == 2
                       X_1 = X_;
                       Y_1 = Y_;                        
                    end                    
                    if R_matched(iter) > R_matched(iter-1) && any(abs(disp_X_Y_med - pixshift) > 5)
                        % If the new matching is better and is
                        % significantly different from the original initial
                        % displacement estimation.
                        X_0 = X_;
                        Y_0 = Y_;
                        pixshift = disp_X_Y_med;
                    else
                        flag_stop = true;
                        X_ = X_0;
                        Y_ = Y_0;
                        rate = R_matched(iter - 1);
                        pixshift = median(X_ - Y_, 1);
                    end
                end
            end
        else
            flag_stop = true;
            pixshift = round(median(disp_X_Y,1));
        end
%% Output
        rate_ = rate;
        nonuniformity = nonuniformity(1 : iter);
    end
%     if isempty(X_) 
%         pixshift = pixshiftinit;
%         return
%     end
end
if isempty(X_)
    pixshift = pixshiftinit;
end
if vis_Q && ~isempty(X_)
    figure;
    scatter3(X_(:,1), X_(:,2), X_(:,3));
    hold on 
    scatter3(Y_(:,1), Y_(:,2), Y_(:,3));
    legend('Tile 1', 'Tile 2');
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
end
end
%% Subfunctions
function [rate,X_,Y_,tY_] = descriptorMatchforz(X,Y,pixshift,iadj,params)
%DESCRIPTORMATCH Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCH(INPUTARGS) Explain usage here
%
% Inputs:
%   X, Y: two 2D real, double marices, specifying the position of the point
%   in the point cloud.
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/23 14:09:29 $	$Revision: 0.1 $
% Copyright: HHMI 2016
% out = [];
% model = params.model;
% debug = params.viz;
%% Initial match based on point drift
[Transform, C] = cpd_register(X,Y,params.opt);
%% check if match is found
% Compute the pairwise euclidean distance between two input array
pD = pdist2(X,Transform.Y);
[aa1,bb1] = min(pD,[],1);
[aa2,bb2] = min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)' < params.projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);
tY_= Transform.Y(keeptheseY(disttrim),:);
% Rate is the ratio of the number of matched pair of distance less than
% projectionThr over the total number of matched pairs
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
% if rate < .5 % dont need to continue
%     [X_,Y_,out] = deal(0);
%     return
% end
Y_ = bsxfun(@minus, Y_, pixshift);
% %%
% % displacement field between follows a field curve on x&y due to
% % optics and deformation curve due to tissue and cut force on z
% dispvec = X_-Y_;
% x = dispvec(:,iadj);
% if iadj==1 % x-neighbor
%     % x : x-displacement
%     % y : y-location
%     y = X_(:,2);
%     bw = [2 220];
% elseif iadj==2 % y-neighbor
%     % x : y-displacement
%     % y : x-location
%     y = X_(:,1);
%     bw=[3 100];
% else % z-neighbor
%     % x : z-displacement
%     % y : y-location (not too much on x as cut is on y direction)
%     y = X_(:,2);
%     bw=[2 220];
% end
% % build a probabilistic model of displacement vectors
% N = 101;
% gridx = linspace(min(x),max(x),N);
% gridy = linspace(min(y),max(y),N);
% [density,bw] = ksdensity2d([x y],gridx,gridy,bw);density=density'/max(density(:));
% [xmin,ix] = min(pdist2(x,gridx'),[],2);
% [ymin,iy] = min(pdist2(y,gridy'),[],2);
% idx = sub2ind([N,N],iy,ix);
% prob_inliers = density(idx)>max(density(idx))*.25;
% x_inline = x(prob_inliers,:);
% y_inline = y(prob_inliers,:);
% %
% % fit curve model
% [~,im] = max(density,[],2);
% sgn = -2*((max(im)==im(end) | max(im)==im(1))-.5);
% pinit = [median(y) sgn*1e-5 median(x)];
% warning off
% out = nlinfit(y_inline, x_inline, model, pinit,optimopts);
% warning on
% % outlier rejection based on parametric model
% xest = feval(model,out,y);
% outliers = abs(x-xest)>2;
% X_ = X_(~outliers,:);
% Y_ = Y_(~outliers,:);
% tY_ = tY_(~outliers,:);
% if debug
%     xgridest = feval(model,out,gridy);
%     figure(100),
%     subplot(2,2,iadj),cla
%     imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
%     axis tight
%     hold on,
%     plot(x,y,'m.')
%     plot(x_inline,y_inline,'mo')
%     plot(x(~outliers),y(~outliers),'gd')
%     plot(xgridest,gridy,'r-')
%     
%     subplot(2,2,4),
%     cla
%     hold on
%     plot3(X(:,1),X(:,2),X(:,3),'b+')
%     plot3(Y(:,1),Y(:,2),Y(:,3),'m.')
%     plot3(Transform.Y(:,1),Transform.Y(:,2),Transform.Y(:,3),'ro')
%     pause(.5)
%     drawnow
%     
% end
end
%% Subfunctions
function [bin_cell_array, varargout] = fun_bin_data_to_idx_list(data)
% fun_bin_data_to_idx_list bin the data according to their values and
% output the corresponding index list
% Input: 
%   data: numerical vector
% Output: 
%   bin_cell_array: cell array, each cell constains a vector, whose
%   components are the indices of the component of data that have the same
%   value. 
%   varargout: unique data value 
% Author: Xiang Ji ( Department of Physics, UC San Diego )
% Nov 28, 2018
if isempty(data)
    bin_cell_array = [];
    varargout{1} = [];
    return;
end
num_data = numel(data);
if ~issorted(data)
    [data, idx_list ]= sort(data, 'ascend');
else
    idx_list = 1 : num_data;
end

bin_size = 0;
bin_data = data(1);
bin_idx = zeros(1, round(num_data/2));
est_num_bin = numel(data);
bin_value_list = zeros(est_num_bin,1);
bin_value_list(1) = data(1);
bin_cell_array = cell(est_num_bin,1);
num_bin = 0;
for idx = 1 : num_data
    tmp_data = data(idx);
    
    if tmp_data == bin_data
        bin_size = bin_size + 1;
        bin_idx(bin_size) = idx_list(idx);
    else
        num_bin = num_bin + 1;
        bin_cell_array{num_bin} = bin_idx(1 : bin_size);
        bin_data = tmp_data;
        bin_value_list(num_bin + 1) = bin_data;
        bin_idx(1) = idx_list(idx);
        bin_size = 1;
    end
end
num_bin = num_bin + 1;
bin_cell_array{num_bin} = bin_idx(1 : bin_size);
bin_cell_array(num_bin + 1 : end) = [];
bin_value_list = bin_value_list(1 : num_bin);
if nargout > 1
    varargout{1} = bin_value_list;
end
end
