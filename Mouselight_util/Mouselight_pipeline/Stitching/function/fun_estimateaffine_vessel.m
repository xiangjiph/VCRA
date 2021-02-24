function [scopeparams,valid_tile_Q] = fun_estimateaffine_vessel(paireddescriptor,neighbors,scopeloc,params,curvemodel,old)
%ESTIMATEAFFINE Summary of this function goes here
%
% [OUTPUTARGS] = ESTIMATEAFFINE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/26 14:29:14 $	$Revision: 0.1 $
% Copyright: HHMI 2017
%%
if nargin<6
    old = 1;
end
beadreg = 0;
checkthese = [1 4 5 7]; % 0 - right - bottom - below
neigs4 = neighbors(:,[1 2 3 4 5]);% left - above - right - bottom
neigs = neighbors(:,checkthese);%[id -x -y +x +y -z +z] format
% beadparams_=beadparams;
dims = params.imagesize;
order = params.order; % power to weight shift
imsize_um = params.imsize_um;

% initialize empty scope params
% beadparamsZmatch_X = [];
% beadparamsZmatch_Y = [];
% beadparamsZmatchdispstage=[];
% TODO: fill below once we get bead images for each scope/objective
% scope = params.scope;
% beadparamsZmatch_X = beadparams.allX{3};
% beadparamsZmatch_Y = beadparams.allY{3};
% if scope==1 % patch
%     beadparams.dispstage{3}=beadparams_.dispstage{3}/150*144.5;
% end
% beadparamsZmatchdispstage = beadparams.dispstage{3};

xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];
% dims is in (x,y,z) order
[sub_x, sub_y] = ndgrid(1:dims(1), 1:dims(2));
xy_1 = [sub_x(:), sub_y(:)];
%%
scopeparams = [];
Ntiles = size(neigs,1);
edges = cell(1,Ntiles);
for iter_tile = 1:Ntiles
    edge = paireddescriptor{iter_tile}.neigs(2:3); % index of 6-neighbor in the positive directions
    counts = paireddescriptor{iter_tile}.count;
    edges{iter_tile} = [[iter_tile;iter_tile],[edge(:)],counts(:)];
end
edges = cat(1,edges{:});
edges(any(isnan(edges),2),:)=[];
G = sparse(edges(:,1),edges(:,2),edges(:,3),Ntiles,Ntiles);
G = max(G,G');
% G is a N-by-N sparese matrix, whose elements are the number of matched
% feature points between tile i and tile j. 
% G is constructed to find the neighboring index of the give tile, since
% neigs only record the neighbors in +x and +y directions. 
%%
parfor_progress(Ntiles)
if old
    skipinds = any(isnan(neigs4(:,[4 5])),2);% neigs4 is NAN in either +X or +Y direction
else
    skipinds = any(isnan(neigs4(:,2:3)),2) & any(isnan(neigs4(:,4:5)),2); % -x/-y is NAN & +x/+y is NAN
end
valid_tile_Q = zeros(1,Ntiles);
%%
try; parfor_progress(0);catch;end
parfor_progress(Ntiles)
for iter_tile = Ntiles : -1 : 1
    %%
    scopeparams(iter_tile).imsize_um = imsize_um; % (x, y, z)
    scopeparams(iter_tile).dims = dims; % (x, y, z)    
    neiginds = find(full(G(iter_tile,:))); 
    theseinds = setdiff(neiginds,paireddescriptor{iter_tile}.neigs(2:3));
    tmp_paired_descriptor = struct;
    % Why X, Y can be non-integer? 
    % Combined all the matched voxel pairs 
    tmp_paired_descriptor.onx.X = round(cat(1, paireddescriptor{iter_tile}.onx.X_skl, paireddescriptor{iter_tile}.onx.X_edge));
    tmp_paired_descriptor.onx.Y = round(cat(1, paireddescriptor{iter_tile}.onx.Y_skl, paireddescriptor{iter_tile}.onx.Y_edge));
    tmp_paired_descriptor.ony.X = round(cat(1, paireddescriptor{iter_tile}.ony.X_skl, paireddescriptor{iter_tile}.ony.X_edge));
    tmp_paired_descriptor.ony.Y = round(cat(1, paireddescriptor{iter_tile}.ony.Y_skl, paireddescriptor{iter_tile}.ony.Y_edge));
    
    if skipinds(iter_tile)
        
    else
        [matched_xyz_1, matched_xyz_2, stage_disp] = deal([]);
        siz = 0; % right/below + left/above
        stgdisp = zeros(3,1); % right/below + left/above
        % right adjacency
        if ~isnan(neigs(iter_tile,2))
            siz(1) = size(tmp_paired_descriptor.onx.X,1);
            matched_xyz_1 = [matched_xyz_1; tmp_paired_descriptor.onx.X];
            matched_xyz_2 = [matched_xyz_2; tmp_paired_descriptor.onx.Y];
            % stage displacement in +x direction
            stgdisp(:,1) = 1000 * ( scopeloc.loc(neigs(iter_tile,2),:) - scopeloc.loc(neigs(iter_tile,1),:) );
            stage_disp = [stage_disp, repelem(1000 * stgdisp(:,1), 1, siz(1))]; 
            valid_tile_Q(iter_tile) = 1;
        end
        if ~isnan(neigs(iter_tile,3))
            siz(end+1) = size(tmp_paired_descriptor.ony.X,1);
            matched_xyz_1 = [matched_xyz_1;tmp_paired_descriptor.ony.X];
            matched_xyz_2 = [matched_xyz_2;tmp_paired_descriptor.ony.Y];
            stgdisp(:,end+1) = 1000 * ( scopeloc.loc(neigs(iter_tile,3),:) - scopeloc.loc(neigs(iter_tile,1),:) );
            stage_disp = [stage_disp, repelem(1000 * stgdisp(:,end), 1, siz(end))];
            valid_tile_Q(iter_tile) = 1;
        end
        % The number of x/y matches can be very unbalanced, since some
        % tiles do not have edge matches
        if isfield(params,'beadparams') && ~isempty(params.beadparams)
            % append bead params if exists
            beadparamsZmatch_X = params.beadparams.allX{3};
            beadparamsZmatch_Y = params.beadparams.allY{3};
            beadparamsZmatchdispstage = params.beadparams.dispstage{3};
        else
            tmp_num = max(1,round(size(matched_xyz_1,1)/2)); % why half? 
            beadparamsZmatch_X = ones(tmp_num,1) * [0 0 250];
            beadparamsZmatch_Y = ones(tmp_num,1) * [0 0 0];
            beadparamsZmatchdispstage = [ones(tmp_num,1) * [0 0 250] * 1e3]';
        end
        siz(end+1) = size(beadparamsZmatch_X,1); % [# matches in x direction, # matches in y direction, an arbitrary number]
        matched_xyz_1 = [matched_xyz_1; beadparamsZmatch_X]; % Add a list of [0, 0, 250] to the matched voxel subscripts in tile 1
        matched_xyz_2 = [matched_xyz_2; beadparamsZmatch_Y]; % Add a list of [0, 0, 250] to the matched voxel subscripts in tile 2
        stage_disp = [stage_disp, beadparamsZmatchdispstage]; % why the stage displacement is 250 when matching the affine transformation in xy direction? 
        matched_xyz_1_FC = matched_xyz_1;
        matched_xyz_2_FC = matched_xyz_2;
        suballX = matched_xyz_1+1;
        suballY = matched_xyz_2+1;
        
        if 1
            % Correct the position of the matched voxels in two adjacent
            % tiles by the curvature model.
            % Q: I don't understand the way with weight is calculated. Why
            % the weight is normalized by p3? 
            [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,iter_tile),order,xy,dims,suballX);
            matched_xyz_1_FC = locs - 1; % coordinate of the matched voxels in tile 1 after correction 
            [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,iter_tile),order,xy,dims,suballY);
            matched_xyz_2_FC = locs - 1; % coordinate of the matched voxels in tile 2 after correction
        else
            cent = squeeze(mean(curvemodel(1:2,1,ineig),3));
            scale = (squeeze(mean(curvemodel(1:2,2,ineig),3)));
            shift = squeeze(mean(curvemodel(1:2,3,ineig),3));
            beta = scale./shift.^order;
            [xshift2D,yshift2D] = shiftxy(xy,cent,beta,order,dims);
            % descriptors are 0 index based
            idxx = sub2ind(dims([2 1]),allX(:,2),allX(:,1));
            xshift = xshift2D(idxx);
            idxy = sub2ind(dims([2 1]),allX(:,2),allX(:,1));
            yshift = yshift2D(idxy);
            allXFC(:,1) = allX(:,1) + xshift;
            allXFC(:,2) = allX(:,2) + yshift;
            idxx = sub2ind(dims([2 1]),allY(:,2),allY(:,1));
            xshift = xshift2D(idxx);
            idxy = sub2ind(dims([2 1]),allY(:,2),allY(:,1));
            yshift = yshift2D(idxy);
            allYFC(:,1) = allY(:,1) + xshift;
            allYFC(:,2) = allY(:,2) + yshift;
        end
        
        Dall = (matched_xyz_2 - matched_xyz_1)';
        DallFC = (matched_xyz_2_FC - matched_xyz_1_FC)';
%         inds = zeros(1,sum(siz(:)));
%         siz12 = [sum(siz(:,1:2),2) siz(:,3)]';
%         idxsub = [0;cumsum(siz12(:))]';
%         for ii=1:2:length(idxsub)-1
%             inds([idxsub(ii)+1:idxsub(ii+1)]) = 1;
%         end
        glS_= stage_disp/Dall;
        glSFC_= stage_disp/DallFC;
        scopeparams(iter_tile).affinegl = glS_;
        scopeparams(iter_tile).affineglFC = glSFC_;
    end
    parfor_progress;
end
parfor_progress(0);
end
