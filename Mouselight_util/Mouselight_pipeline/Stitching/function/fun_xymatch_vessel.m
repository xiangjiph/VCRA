function [paireddescriptor,R,curvemodel] = fun_xymatch_vessel(descriptors,neigs,scopeloc,params)
%ESTIMATESCOPEPARAMETERS Summary of this function goes here
%
% [OUTPUTARGS] = ESTIMATESCOPEPARAMETERS(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/12 10:38:28 $	$Revision: 0.1 $
% Copyright: HHMI 2016
%%
% addpath(genpath('./thirdparty'))
debug_mode = false;
% res = 0; 
% viz=0;
fignum = 101;
projectionThr = 20; % distance between target and projected point has to be less than this number
tile_size_xyz = params.imagesize;
imsize_um = params.imsize_um;
min_num_skl_match = 50;
min_num_edge_match = 300;
% slid = [[75 960];[0 dims(2)];[0 dims(3)]];
% expensionshift = [0 0 20]; % HEURISTICS:: tissue expends, so overlap is bigger between tiles
%% Set parameters for feature matching
model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
optimopts = statset('nlinfit');
optimopts.RobustWgtFun = 'bisquare';
% opt.method='nonrigid_lowrank';
opt.method='nonrigid';
opt.beta=6;            % the width of Gaussian kernel (smoothness), higher numbers make transformation more stiff
opt.lambda=16;          % regularization weight
opt.viz=0;              % show every iteration
opt.outliers=0.9;       % use 0.7 noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
%     opt.max_it=100;         % max number of iterations
%     opt.tol=1e-10;          % tolerance
matchparams.model = model;
matchparams.optimopts = optimopts;
matchparams.opt = opt;
matchparams.projectionThr = projectionThr;
matchparams.debug = debug_mode;
matchparams.viz = debug_mode;
matchparams.fignum = fignum;
matchparams.opt.beta=2;
% matchparams.opt.method = 'nonrigid';
% What does the following values mean? 
matchparams.init(1,:)=[765 1e-5 867];
matchparams.init(2,:)=[537 -1e-5 1445];
%% Initialization
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% indicies are 1 based,e.g. x = 1:dims(1), not 0:dims(1)-1
% xyz_umperpix = zeros(size(neigs,1),3);
curvemodel = zeros(3,3,size(neigs,1));
R = zeros(3,3,size(neigs,1));

paireddesctemp = [];
paireddesctemp{1}.valid = false;
paireddesctemp{1}.pixshift_mask_fft = [];
paireddesctemp{1}.match_rate_mask_fft = [];
paireddesctemp{1}.X_skl = [];
paireddesctemp{1}.Y_skl = [];
paireddesctemp{1}.count_skl = [];
paireddesctemp{1}.pixshift_skl = [];
paireddesctemp{1}.match_rate_skl = [];
paireddesctemp{1}.X_edge = [];
paireddesctemp{1}.Y_edge = [];
paireddesctemp{1}.pixshift_edge = [];
paireddesctemp{1}.match_rate_edge = [];
paireddesctemp{1}.count_edge = [];
paireddesctemp{1}.count = [];
paireddesctemp{1}.pixshift_stage = [];
paireddesctemp{1}.fit_inliers_ratio = 0;
paireddesctemp{2} = paireddesctemp{1};


paireddescriptor = cell(size(neigs,1),1);
for ix = 1:size(neigs,1)
    paireddescriptor{ix}.neigs = neigs(ix,:);
    paireddescriptor{ix}.ony = paireddesctemp{1};
    paireddescriptor{ix}.onx = paireddesctemp{1};
    paireddescriptor{ix}.count = [];
end

Ntiles = size(neigs,1);
%%
try parfor_progress(0);catch;end
parfor_progress(Ntiles)
num_core = feature('numcores');
num_process = round(num_core * 0.75);
use_existing_xy_match = true;
update_outlier_rejectionQ = false;
parfor (ineig = 1 : Ntiles, num_process)
    %% load descriptor pairs X (center) - Y (adjacent tile)
    fprintf('Processing %d / %d\n', ineig, Ntiles);
    idx_center = neigs(ineig,1);
    des_center = descriptors{idx_center};
    if isempty(des_center);continue;end
    mout = zeros(3,3);
    paireddescriptor_ = paireddesctemp;
    R_ = zeros(3, 3); 

    for iadj = 1 : size(neigs,2)-2 %1:x-overlap, 2:y-overlap, 3:z-overlap
%% Check if the tile has been processed or not
        read_precomputed_matchesQ = false;
        if use_existing_xy_match && ~isempty(des_center.record) && isfield(des_center.record, 'fp_descriptor')
            tmpfp = strrep(des_center.record.fp_descriptor, 'stage_2_descriptor_output', ...
                'stage_3_point_match_output');
            % For local test
%             tmpfp = strrep(tmpfp, '/nrs/mouselight/pipeline_output/', '/data/Vessel/ML_stitching/');
            
            tmp_folder = fileparts(tmpfp);
            switch iadj
                case 1
                    tmpfp = fullfile(tmp_folder, 'match-X.mat');
                case 2
                    tmpfp = fullfile(tmp_folder, 'match-Y.mat');
                case 3
%                     tmpfp = fullfile(tmp_folder, 'match-Z.mat');
            end
            
            if isfile(tmpfp)
                tmp_str = load(tmpfp);
                read_precomputed_matchesQ = true;
                tmp_str = tmp_str.tmp_xy_match_result;
                paireddescriptor_{iadj} = tmp_str.paireddescriptor;
                R_(:,iadj) = tmp_str.R_;
                mout(iadj,:) = tmp_str.mout
                
                if ~update_outlier_rejectionQ
                    continue;
                end                
            end
        end
%% Estimate stage shift
        % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
        idx_adj =  neigs(ineig,iadj+1);
        if isnan(idx_adj);continue;end
        stgshift = 1000*(scopeloc.loc(idx_adj,:)-scopeloc.loc(idx_center,:));
        pixshift_stage = round(stgshift.*(tile_size_xyz-1)./(imsize_um));
        paireddescriptor_{iadj}.pixshift_stage = pixshift_stage;
        des_adj = descriptors{idx_adj};
        if isempty(des_adj);continue;end
        
        if ~read_precomputed_matchesQ 
%% Intensity based masked fft registration 
% We never use this information computed here. 
            [paireddescriptor_{iadj}.pixshift_mask_fft, paireddescriptor_{iadj}.match_rate_mask_fft] = deal([]);
%         tic
%         [paireddescriptor_{iadj}.pixshift_mask_fft, paireddescriptor_{iadj}.match_rate_mask_fft] = ...
%             fun_masked_fft_match_vessel(des_center, des_adj, pixshift_stage, iadj, debug_mode);
%         toc
%% Point cloud registration
            if ~isempty(des_center.skl_sub) && ~isempty(des_adj.skl_sub)
                des_center_skl = cat(2, util.correctTiles(des_center.skl_sub), des_center.skl_label(:));
                des_adj_skl = cat(2, util.correctTiles(des_adj.skl_sub), des_adj.skl_label(:));
                matchparams_skl = struct;
                % Parameter for CPD
                matchparams_skl.opt.method='nonrigid';
                matchparams_skl.opt.beta = 6;            % the width of Gaussian kernel (smoothness), higher numbers make transformation more stiff
                matchparams_skl.opt.lambda = 16;          % regularization weight
                matchparams_skl.opt.viz = 0;              % show every iteration
                matchparams_skl.opt.outliers = 0.9;       % use 0.9 noise weight
                matchparams_skl.opt.fgt = 0;              % do not use FGT (default)
                matchparams_skl.opt.normalize = 1;        % normalize to unit variance and zero mean before registering (default)
                matchparams_skl.opt.corresp = 1;          % compute correspondence vector at the end of registration (not being estimated by default)
                % Parameters for skeleton matching
                matchparams_skl.max_num_desc = 1e4;     % Maximum number of skeleton voxels for matching
                matchparams_skl.scan_pixshift_Q = false; % Do not search in the x-y direction, as the stage should not be that off. 
                %         matchparams.projectionThr = projectionThr;
                matchparams_skl.projectionThr = 5;
                matchparams_skl.viz = debug_mode;

                [X_skel,Y_skel,rate_skl, pixshift_skl,~] = searchpair_vessel(des_center_skl,des_adj_skl,...
                    pixshift_stage, iadj, tile_size_xyz, matchparams_skl);
    %           Save information: no matter how badly the matches is, saving
    %           the matching faciliate the subsequent outlier rejection test. 
                paireddescriptor_{iadj}.X_skl = util.correctTiles(X_skel,tile_size_xyz);
                paireddescriptor_{iadj}.Y_skl = util.correctTiles(Y_skel,tile_size_xyz);
                paireddescriptor_{iadj}.count_skl = size(paireddescriptor_{iadj}.X_skl, 1);
                paireddescriptor_{iadj}.pixshift_skl = pixshift_skl;
                paireddescriptor_{iadj}.match_rate_skl = rate_skl;
            else
                [X_skel, Y_skel, rate_skl, pixshift_skl ] = deal([]);
            end        
%% Edge match
            if isfield(des_center.record, 'compute_edge') &&  des_center.record.compute_edge && ...
                    isfield(des_adj.record, 'compute_edge') && des_adj.record.compute_edge
                disp('Vessel edge CPD');
                desc1_edge = cat(2, des_center.edge_sub, des_center.edge_gradient);
                desc2_edge = cat(2, des_adj.edge_sub, des_adj.edge_gradient);
                desc1_edge = correctTiles(desc1_edge, tile_size_xyz);
                desc2_edge = correctTiles(desc2_edge, tile_size_xyz);
    %             tic
                try
                    [X_edge, Y_edge, rate_edge, pixshift_edge] = fun_searchpair_vessel_edges(desc1_edge, desc2_edge, pixshift_stage);
                catch ME
                    if (strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded'))
                        % Probably empty tile.
                        X_edge = [];
                        Y_edge = [];
                        rate_edge = [];
                        pixshift_edge = [];
                    elseif (strcmp(ME.identifier, 'MATLAB:nomen'))
                        warning('Out of memory. Set X_edge, etc to [].');
                        % Out of memory error. Should not happen for normal
                        % tile...since the edge voxels have readly been merged
                        % before pwdist2
                        X_edge = [];
                        Y_edge = [];
                        rate_edge = [];
                        pixshift_edge = [];
                    else
                        rethrow(ME);
                    end
                end
    %             toc
                % Save information
                paireddescriptor_{iadj}.X_edge = util.correctTiles(X_edge, tile_size_xyz);
                paireddescriptor_{iadj}.Y_edge = util.correctTiles(Y_edge, tile_size_xyz);
                paireddescriptor_{iadj}.pixshift_edge = pixshift_edge;
                paireddescriptor_{iadj}.match_rate_edge = rate_edge;
                paireddescriptor_{iadj}.count_edge = size(paireddescriptor_{iadj}.X_edge, 1);
            else
                [X_edge, Y_edge, rate_edge, pixshift_edge] = deal([]);
            end
        else
            % Use the pre-computed matches
            X_skel = util.correctTiles(paireddescriptor_{iadj}.X_skl, tile_size_xyz);
            Y_skel = util.correctTiles(paireddescriptor_{iadj}.Y_skl, tile_size_xyz);
            rate_skl = paireddescriptor_{iadj}.match_rate_skl;
            
            X_edge = util.correctTiles(paireddescriptor_{iadj}.X_edge, tile_size_xyz);
            Y_edge = util.correctTiles(paireddescriptor_{iadj}.Y_edge, tile_size_xyz);
            rate_edge = paireddescriptor_{iadj}.match_rate_edge;
        end
%% Outlier rejection 
        if size(X_skel, 1) < min_num_skl_match || rate_skl < 0.8
            X_skel = [];
            Y_skel = [];
        end
        if size(X_edge, 1) < min_num_edge_match || rate_edge < 0.8
            X_edge = [];
            Y_skel = [];
        end
%% Get field curvature model
% For iadj = 1, i.e. x-neighbor tile, Use model = @(p,y) p(3) -
% p(2).*((y-p(1)).^2) to fit the pairwise displacement in x-direction
% between X_skel and Y_skel with respect to the y-coordinate of X_skel. The
% fitted model is used to reject outliers in the input voxel pairs. If the
% percentage of the outlier is less than 25%, set valid = 1, and the
% outliers are removed and the resulting voxel pairs are output to X_ and
% Y_ respectively. out is the fitting parameters.
%         [X_,Y_,out,valid] = match.fcestimate(X_skel, Y_skel, iadj, matchparams);
%         [X_edge_fc,Y_edge_fc,out_edge,valid_edge] = match.fcestimate(X_edge, Y_edge, iadj, matchparams);
% The fitting parameters from skeleton and edge are not very consistent 
        % Merge two matched dataset and randomly uniformly sample: 
        if ~isempty(X_skel) || ~isempty(X_edge)
            tmp_grid_size = 1;
            sub_sample_Q = true;
            tmp_merge_dim = setdiff([1,2], iadj);
            X_matched = cat(1, X_skel, X_edge);
            Y_matched = cat(1, Y_skel, Y_edge);
            if sub_sample_Q
                sub_min = min(X_matched(:, tmp_merge_dim), [], 1);
    %             sub_max = max(X_matched(:, tmp_merge_dim), [], 1);
    %             tmp_image_size = sub_max - sub_min + 1;
                tmp_bbox_sub = ceil((1 + bsxfun(@minus, X_matched(:, tmp_merge_dim), sub_min))./ tmp_grid_size);
                [tmp_bin_cell_array, tmp_bin_ind] = fun_bin_data_to_idx_list(tmp_bbox_sub);
                num_sample = numel(tmp_bin_ind);
                sample_ind = zeros(num_sample, 1);
                for iter_bin = 1 : numel(tmp_bin_ind)
                    tmp_ind_list = tmp_bin_cell_array{iter_bin};
                    if length(tmp_ind_list) > 1
                        sample_ind(iter_bin) = randsample(tmp_ind_list,1);
                    else
                        sample_ind(iter_bin) = tmp_ind_list;
                    end
                end
                X_sampled = X_matched(sample_ind, :);
                Y_sampled = Y_matched(sample_ind, :);
            end

            fc_params = fun_fcestimate_vessel(X_sampled, Y_sampled, iadj, matchparams);
            if fc_params.valid
                % Use the field curvature fitting to reject outliers in the
                % edge and skeleton matches.
                fc_params.viz = false;
                [X_edge, Y_edge] = fun_matched_pair_rejection_by_fcestimate_vessel(X_edge, Y_edge, fc_params);
%                  The interpretation of the fitting parameters is not very
%                  clear.
                [X_skel, Y_skel] = fun_matched_pair_rejection_by_fcestimate_vessel(X_skel, Y_skel, fc_params);
            end
            mout(iadj,:) = fc_params.fitting_params;
            paireddescriptor_{iadj}.valid = fc_params.valid;
            paireddescriptor_{iadj}.fit_inliers_ratio = fc_params.inliers_ratio; 
        end
%% Select and merge matched voxels
% To avoid the unbalanced contribution from the edge skeleton, further
% sample the matched edge skeleton? 
% Trust the skeleton more than the edge - Actually the problem is not this
% tile, since both the skeleton and the edge are available. 
% Merge the list and sample evenly across the occupied space 
% Sample the edge and skeleton saperately
        sample_block_scale = 15;
        sample_block_size = [3, 3, 1] .* sample_block_scale;
        num_sample_per_block = ceil(sample_block_scale/3);
        matched_1_sub = cat(1, X_skel, X_edge);
        matched_2_sub = cat(1, Y_skel, Y_edge);    
        [matched_1_sub, sampled_ind] = fun_uniform_sample_points_in_space(matched_1_sub, sample_block_size, num_sample_per_block, 'random');
        matched_2_sub = matched_2_sub(sampled_ind, :);
%         downsample_edgeQ = true;
%         downsample_skelQ = false;
%         if ~isempty(X_skel) && downsample_skelQ
%             [X_skel_sampled, sampled_ind] = fun_uniform_sample_points_in_space(X_skel, sample_block_size, num_sample_per_block, 'random');
%             Y_skel_sampled = Y_skel(sampled_ind, :);
%         else
%             X_skel_sampled = X_skel;
%             Y_skel_sampled = Y_skel;
%         end
%         if ~isempty(X_edge) && downsample_edgeQ
%             [X_edge_sampled, sampled_ind] = fun_uniform_sample_points_in_space(X_edge, sample_block_size, num_sample_per_block, 'random');
%             Y_edge_sampled = Y_edge(sampled_ind, :);
%         else
%             X_edge_sampled = X_edge;
%             Y_edge_sampled = Y_edge;
%         end
%         matched_1_sub = cat(1, X_skel_sampled, X_edge_sampled);
%         matched_2_sub = cat(1, Y_skel_sampled, Y_edge_sampled);            
%% Save information
        paireddescriptor_{iadj}.count_before_sample = paireddescriptor_{iadj}.count_edge + paireddescriptor_{iadj}.count_skl;
        paireddescriptor_{iadj}.sample_block_size = sample_block_size;
        
        paireddescriptor_{iadj}.X = util.correctTiles(matched_1_sub,tile_size_xyz);
        paireddescriptor_{iadj}.Y = util.correctTiles(matched_2_sub,tile_size_xyz);
        paireddescriptor_{iadj}.count = size(matched_1_sub,1);        
        %R(:,iadj,ineig) = round(median(X_-Y_));
        % !!! Here X_ and Y_ have been flipped, so the R_ is exactly the
        % opposite of (X_ - Y_) before flip
        % For now, choose the skeleton voxel 
        R_(:,iadj) = round(median(X_skel - Y_skel));
        % Save partial information
        if ~isempty(des_center.record) && isfield(des_center.record, 'fp_descriptor')
            tmp_xy_match_result = struct;            
            tmp_xy_match_result.R_ = R_(:,iadj);
            tmp_xy_match_result.mout = mout(iadj,:);
            tmp_xy_match_result.paireddescriptor = paireddescriptor_{iadj};
            
            tmpfp = strrep(des_center.record.fp_descriptor, 'stage_2_descriptor_output', ...
                'stage_3_point_match_output');
            tmp_folder = fileparts(tmpfp);
            if iadj == 1
               tmpfp = fullfile(tmp_folder, 'match-X.mat');
            elseif iadj == 2
                tmpfp = fullfile(tmp_folder, 'match-Y.mat');
            elseif iadj == 3
%                 tmpfp = fullfile(tmp_folder, 'match-Z.mat');
            end
            if ~isfolder(tmp_folder)
                warning('Target folder does not exist');
                mkdir(tmp_folder);
                unix(sprintf('chmod g+rw %s',tmp_folder));
            end
            parsave(tmpfp, tmp_xy_match_result);
            unix(sprintf('chmod g+rw %s',tmpfp));
        end
    end
    R(:,:,ineig) = R_;
    curvemodel(:,:,ineig) = mout;
    paireddescriptor{ineig}.onx = paireddescriptor_{1};
    paireddescriptor{ineig}.ony = paireddescriptor_{2};
    paireddescriptor{ineig}.count = paireddescriptor_{1}.count + paireddescriptor_{2}.count;
    parfor_progress;
end
parfor_progress(0);
end
