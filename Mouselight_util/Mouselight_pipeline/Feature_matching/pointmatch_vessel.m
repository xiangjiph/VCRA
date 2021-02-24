function varargout = pointmatch_vessel(tile1,tile2,acqusitionfolder1,acqusitionfolder2,outfold,pixshift,ch,maxnumofdesc,exitcode)
% pointmatch_vessel finds the corresponding position of the same object in
% the adjacent tiles by the following steps:
% 1. Use the microscope stage position information to estimate the relative
% displacement between two tiles. Load the descriptors from
% vesselDescriptor.
% 2. If two tiles both contain large vessels near the boundary and the two
% tiles are adjacent in z direction, use masked intensity based FFT
% registration to estiamte the translation transformation between two
% tiles.
% ( Need to decide when to use the displacement from masked FFT to update
% the pixshift for the following computation)
% 3. Apply Coherent Point Drift algorithm to find the correspondence
% between the the skeleton of the vessel from two tiles. Function used:
% searchpair_vessel
% 4. If both tiles contain large vessels near the boundary, use the vessel
% edge for point cloud registration. Function used: fun_searchpair_vessel_edges
% 5. Record and save all the registration result.
%
% Modified from Erhan Bas's pointmatch by Xiang Ji (xiangji.ucsd@gmail.com)
% Date: Dec 14, 2018
% dbstop if error
varargout{1} = -1;
%% Complied file setting
compiledfunc = '/groups/mousebrainmicro/home/jix/Documents/Github/compiledfunctions/pointmatch_vessel/pointmatch_vessel';
if ~exist(fileparts(compiledfunc),'dir')
    mkdir(fileparts(compiledfunc));
    mfilename_ = mfilename('fullpath');
    % unix(sprintf('mcc -m -v -R -singleCompThread %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions')))
    complie_cmd_str = sprintf('mcc -m -R -nodisplay -v %s -d %s -a %s',mfilename_,fileparts(compiledfunc),fullfile(fileparts(mfilename_),'functions'));
    unix(complie_cmd_str)
    unix(sprintf('chmod g+x %s',compiledfunc))
    return;
end
%%
disp('Start processing');
debug_mode = false;
overwrite_existing_matching = false;
%% Path setting
if ~isdeployed
    addpath(genpath('./functions'))
end
%% Default input
if nargin<1
    rawfolder = '/groups/mousebrainmicro/mousebrainmicro/data/';
    classifierfolder = '/nrs/mouselight/pipeline_output/';
    sample = '2018-08-15-skeltest';
    tileid1 = '/2018-08-18/00/00152';
    tileid2 = '/2018-08-18/00/00442';
    tile1 = fullfile(classifierfolder,sample,'/stage_2_descriptor_output',tileid1);
    tile2 = fullfile(classifierfolder,sample,'/stage_2_descriptor_output',tileid2);
    acqusitionfolder1 = fullfile(rawfolder,sample,'Tiling',tileid1);
    acqusitionfolder2 = fullfile(rawfolder,sample,'Tiling',tileid2);
end

if nargin < 5
    outfold = tile1;
    pixshift = '[0 0 0]';
    ch='1';
    maxnumofdesc = 1e4;
    exitcode = 0;
elseif nargin < 6
    pixshift = '[0 0 0]';
    ch='1';
    maxnumofdesc = 1e4;
    exitcode = 0;
elseif nargin < 7
    ch='1';
    maxnumofdesc = 1e4;
    exitcode = 0;
elseif nargin < 8
    maxnumofdesc = 1e4;
    exitcode = 0;
elseif nargin < 9
    exitcode = 0;
end
%% Variable conversion
if ischar(pixshift)
    pixshift = eval(pixshift); % pass initialization
end
if ischar(maxnumofdesc)
    maxnumofdesc = str2double(maxnumofdesc);
end
if ischar(exitcode)
    exitcode = str2double(exitcode);
end
%% Use the scope position to estimate the pixel shift
if length(ch)>1
    ch_desc={ch(1),ch(2)};
else
    ch_desc={ch};
end
varargout{1} = exitcode;
%%
tile_size_xyz = [1024,1536,251];
projectionThr = 5;

tag = 'XYZ';
scopefile1 = readScopeFile(acqusitionfolder1);
scopefile2 = readScopeFile(acqusitionfolder2);
imsize_um = [scopefile1.x_size_um,scopefile1.y_size_um,scopefile1.z_size_um];
% estimate translation
gridshift = ([scopefile2.x scopefile2.y scopefile2.z]-[scopefile1.x scopefile1.y scopefile1.z]);
iadj = find(gridshift);
% Stage shift in micron
stgshift = 1000*([scopefile2.x_mm scopefile2.y_mm scopefile2.z_mm]-[scopefile1.x_mm scopefile1.y_mm scopefile1.z_mm]);
% Conver the stage shift to voxel distance
if all(pixshift==0)
    pixshift = round(stgshift .* (tile_size_xyz-1) ./imsize_um);
end
if pixshift(3) >= 235
    warning('Estimated pixel shift in z direction is greater than 249');
    % Apply manual correction 
    pixshift(3) = 235;
end
pixshift_0_edge = round(pixshift);
paireddescriptor = struct;
paireddescriptor.pixshift_stage = pixshift;
[paireddescriptor.exist_blv, paireddescriptor.matchrate_edge, ...
    paireddescriptor.X_edge, paireddescriptor.Y_edge, ...
    paireddescriptor.X_edge_int, paireddescriptor.Y_edge_int, ...
    paireddescriptor.pixshift_edge, paireddescriptor.pixshift_mask_fft, ...
    paireddescriptor.matchrate_mask_fft, paireddescriptor.pixshift_skl,...
    paireddescriptor.matchrate,paireddescriptor.X_skl, paireddescriptor.Y_skl, ...
    paireddescriptor.X_skl_int, paireddescriptor.Y_skl_int, ...
    paireddescriptor.uni] = deal([]);
%%
% check if input exists
% For the compatibility of the original pipeline
[~, ~, tmp_ext] = fileparts(tile1);
if isempty(tmp_ext) %~strcmp(tmp_ext, '.mat')
    tmp = dir(fullfile(tile1, '*desc.0.mat'));
    tile1 = fullfile(tmp.folder, tmp.name);
    tmp = dir(fullfile(tile2, '*desc.0.mat'));
    tile2 = fullfile(tmp.folder, tmp.name);
end

if ~isfile(tile1) || ~isfile(tile2)
    disp('Data for at least 1 tile does not exist. Skip');
    rate_ = 0;
    X_skel = [];
    Y_skel = [];
    uni = 0;
else
    descriptor_1 = load(tile1);
    descriptor_2 = load(tile2);
    disp('Finish loading descriptors');
    %% for debug - fixing the bug in vessel descriptor, should be removed later
    %     if isempty(descriptor_1.record)
    %         disp('Missing field: fp_image. Infer image path automatically');
    %         descriptor_1.record.fp_image = strrep(strrep(tile1, 'stage_2_descriptor_output', 'raw_data'), 'descriptor.mat', '.tif');
    %     end
    %     if isempty(descriptor_2.record)
    %         disp('Missing field: fp_image. Infer image path automatically');
    %         descriptor_2.record.fp_image = strrep(strrep(tile2, 'stage_2_descriptor_output', 'raw_data'), 'descriptor.mat', '.tif');
    %     end
    %% Intensity based masekd fft registration
    % im_1 = deployedtiffread(descriptor_1.record.fp_image);
    % im_1 = flip(flip(im_1, 1), 2);
    % im_2 = deployedtiffread(descriptor_2.record.fp_image);
    % im_2 = flip(flip(im_2, 1), 2);
    % implay(im_1);
    % implay(im_2);
    %% Skeleton point registration
    if isfield(descriptor_1, 'skl_sub') && ~isempty(descriptor_1.skl_sub) && isfield(descriptor_2, 'skl_sub') && ~isempty(descriptor_2.skl_sub)
        desc1_skel = cat(2, correctTiles(descriptor_1.skl_sub,tile_size_xyz), descriptor_1.skl_label(:));
        desc2_skel = cat(2, correctTiles(descriptor_2.skl_sub,tile_size_xyz), descriptor_2.skl_label(:));
        matchparams = modelParams(projectionThr, debug_mode); % Setting parameters for the matching algorithm
        matchparams.max_num_desc = maxnumofdesc;
        matchparams.scan_z_shift_Q = true;
        matchparams.vis = false;
        matchparams.scan_pixshift_Q = true;
        matchparams.remove_descriptor_by_pdist2 = true;
        matchparams.pixshift_search_shift_z = [15, 30, 50];
        matchparams.min_pts_good_match = 100;
        matchparams.min_pts_no_mask_fft = 500;
        if length(iadj)~=1 || max(iadj)>3
            error('not 6 direction neighbor')
        end
        %% MATCHING
        disp('Vessel skeleton CPD');
        tic
        [X_skel, Y_skel, rate_, pixshift_skl, nonuniformity] = searchpair_vessel(desc1_skel, desc2_skel, pixshift, iadj, tile_size_xyz, matchparams);
        toc
        if size(X_skel, 1) < matchparams.min_pts_no_mask_fft% If the matching is not good enough, use the mask fft to initialize the position and match again
            %     disp('Masked FFT translation registration');
            disp('Try to refine initial position estimation by masked fft translational registration');
            tic
            [paireddescriptor.pixshift_mask_fft, paireddescriptor.matchrate_mask_fft] = fun_masked_fft_match_vessel(descriptor_1, descriptor_2, pixshift, iadj, debug_mode);
            toc
            if ~isempty(paireddescriptor.matchrate_mask_fft) && ...
                    (any(abs(paireddescriptor.pixshift_mask_fft - pixshift) > [10, 10, 30]) || ...
                    any(abs(pixshift_skl - paireddescriptor.pixshift_mask_fft) > 10))
                % If the displacement found by mask fft is far from the
                % stage estimation
                disp('Recompute the point cloud registration');
                [X_skel_1, Y_skel_1, rate_1, pixshift_skl_1, nonuniformity_1] = ...
                    searchpair_vessel(desc1_skel, desc2_skel, paireddescriptor.pixshift_mask_fft, iadj, tile_size_xyz, matchparams);
                if size(X_skel_1, 1) > size(X_skel, 1)
                    % If more matching is found, use the new one. 
                    X_skel = X_skel_1;
                    Y_skel = Y_skel_1;
                    rate_ = rate_1;
                    pixshift_skl = pixshift_skl_1;
                    nonuniformity = nonuniformity_1;
                end
            end
        end
        if ~isempty(X_skel)
            [~, ~, matched_idx_1] = intersect(X_skel, desc1_skel(:, 1:3), 'rows', 'stable');
            if numel(matched_idx_1) ~= size(X_skel, 1)
                warning('Not all the points in X_skel are in desc1_skel');
            end
            paireddescriptor.X_skl_int = descriptor_1.skl_int(matched_idx_1);
            
            [~, ~, matched_idx_2] = intersect(Y_skel, desc2_skel(:, 1:3), 'rows', 'stable');
            if numel(matched_idx_2) ~= size(Y_skel, 1)
                warning('Not all the points in Y_skel are in desc2_skel');
            end
            paireddescriptor.Y_skl_int = descriptor_2.skl_int(matched_idx_2);
            X_skel = correctTiles(X_skel,tile_size_xyz);
            Y_skel = correctTiles(Y_skel,tile_size_xyz);
        end
        uni = mean(nonuniformity)<=.5;
        paireddescriptor.pixshift_skl = pixshift_skl;
        paireddescriptor.matchrate = rate_;
        paireddescriptor.X_skl = X_skel;
        paireddescriptor.Y_skl = Y_skel;
        paireddescriptor.uni = uni;
        if rate_ > 0.8 && size(X_skel, 1) > matchparams.min_pts_good_match
            pixshift_0_edge = pixshift_skl;
        end
    else
        [X_skel, Y_skel] = deal([]);
    end
    %% If both the descriptor contains boundary large vessels
    
    if isfield(descriptor_1.record, 'compute_edge') && descriptor_1.record.compute_edge && ...
            ~isempty(descriptor_1.edge_sub) && isfield(descriptor_2.record, 'compute_edge') &&...
            descriptor_2.record.compute_edge && ~isempty(descriptor_2.edge_sub)
        desc1_edge = descriptor_1.edge_sub;
        desc2_edge = descriptor_2.edge_sub;
        desc1_edge = cat(2, correctTiles(desc1_edge, tile_size_xyz), descriptor_1.edge_gradient);
        desc2_edge = cat(2, correctTiles(desc2_edge, tile_size_xyz), descriptor_2.edge_gradient);
        disp('Vessel edge CPD');
        %         tic
        try
            [X_edge, Y_edge, rate_edge, pixshift_edge] = fun_searchpair_vessel_edges(desc1_edge, desc2_edge, pixshift_0_edge);
        catch ME
            if (strcmp(ME.identifier, 'MATLAB:array:SizeLimitExceeded'))
                % Probably empty tile.
                X_edge = [];
                Y_edge = [];
                rate_edge = [];
                pixshift_edge = [];
            elseif (strcmp(ME.identifier, 'MATLAB:nomen'))
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
        %         toc
        if ~isempty(X_edge)
            [~, ~, matched_idx_1] = intersect(round(X_edge), desc1_edge(:, 1:3), 'rows', 'stable');
            paireddescriptor.X_edge_int = descriptor_1.edge_int(matched_idx_1);
            [~, ~, matched_idx_2] = intersect(round(Y_edge), desc2_edge(:, 1:3), 'rows', 'stable');
            paireddescriptor.Y_edge_int = descriptor_2.edge_int(matched_idx_2);
            
            X_edge = correctTiles(X_edge, tile_size_xyz);
            Y_edge = correctTiles(Y_edge, tile_size_xyz);
        end
        paireddescriptor.exist_blv = true;
        paireddescriptor.matchrate_edge = rate_edge;
        paireddescriptor.X_edge = X_edge;
        paireddescriptor.Y_edge = Y_edge;
        paireddescriptor.pixshift_edge = pixshift_edge;
    else
        paireddescriptor.exist_blv = false;
    end
end
%% Downsampling
paireddescriptor.X = cat(1, paireddescriptor.X_skl, paireddescriptor.X_edge);
paireddescriptor.Y = cat(1, paireddescriptor.Y_skl, paireddescriptor.Y_edge);
% Treat this output as raw data. Do the downsampling later, do not throw
% out any information at this stage.
%% Output
if nargout == 2
    varargout{2} = paireddescriptor;
end
if nargin>4
    if ~isfolder(outfold)
        %         warning('Output folder does not exist. Create folder');
        mkdir(outfold);
        unix(sprintf('chmod g+rw %s',outfold));
    end
    if ~isempty(X_skel)
        %x:R, y:G, z:B
        col = median(Y_skel - X_skel, 1)+128;
        col = max(min(col,255),0);
        outpng = zeros(105,89,3);
        outpng(:,:,1) = col(1);
        outpng(:,:,2) = col(2);
        outpng(:,:,3) = col(3);
        if exist(fullfile(outfold,'Thumbs.png'),'file')
            % -f to prevent prompts
            unix(sprintf('rm -f %s',fullfile(outfold,'Thumbs.png')));
        end
        imwrite(outpng,fullfile(outfold,'Thumbs.png'))
        unix(sprintf('chmod g+rw %s',fullfile(outfold,'Thumbs.png')));
    end
    %%
    if isempty(outfold)
        varargout{2} = paireddescriptor;
    else
        % if isempty(rate_); val=0;elseif rate_<1;val=0;else;val=1;end
        outputfile = fullfile(outfold,sprintf('match-%s.mat',tag(iadj))); % append 1 if match found
        
        disp('Write matching result to folder');
        %check if file exist
        if exist(outputfile,'file') && ~overwrite_existing_matching
            % if main match exists, crete a versioned one
            outputfile1 = fullfile(outfold,sprintf('match-%s-1.mat',tag(iadj))); % append 1 if match found
            save(outputfile1,'paireddescriptor','scopefile1','scopefile2')
            unix(sprintf('chmod g+rw %s',outputfile1));
        else
            save(outputfile,'paireddescriptor','scopefile1','scopefile2')
            unix(sprintf('chmod g+rw %s',outputfile));
        end
    end
end

end
% Subfunctions
function [Iout] = deployedtiffread(fileName,slices)
%DEPLOYEDTIFFREAD Summary of this function goes here
%
% [OUTPUTARGS] = DEPLOYEDTIFFREAD(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/21 12:26:16 $	$Revision: 0.1 $
% Copyright: HHMI 2015
warning off
info = imfinfo(fileName, 'tif');
if nargin<2
    slices = 1:length(info);
end
wIm=info(1).Width;
hIm=info(1).Height;
numIm = numel(slices);
Iout  = zeros(hIm, wIm, numIm,'uint16');

for i=1:numIm
    Iout(:,:,i) = imread(fileName,'Index',slices(i),'Info',info);
end

end

function output = crop_bbox3(data, bbox_parameters, bbox_order)
% CROP_BBOX3 crops part of the array DATA according to the given bounding
% box parameters.
% default bbox_parameters = [ul1, ul2, ul3, l1, l2, l3]
% matlab's regionpros3 output bbox = [ul2, ul1, ul3, l2, l1, l3]
if nargin < 3
    bbox_order = 'default';
    warning('bbox_parameters order not specify. Option: default/ regionprop');
end
if ~iscell(bbox_order)
    bbox_parameters = num2cell(round(bbox_parameters));
end
switch bbox_order
    case {'default'}
        if length(bbox_parameters) <=4
            [ul1, ul2, ul3, l1] = bbox_parameters{:};
            l2 = l1;
            l3 = l1;
        else
            [ul1, ul2, ul3, l1, l2, l3] = bbox_parameters{:};
        end
    case {'regionprop'}
        if length(bbox_parameters) <=4
            [ul2, ul1, ul3, l1] = bbox_parameters{:};
            l2 = l1;
            l3 = l1;
        else
            [ul2, ul1, ul3, l2, l1, l3] = bbox_parameters{:};
        end
end
output = data(ul1:ul1+l1-1, ul2:ul2+l2-1, ul3:ul3+l3-1);
end

function desc = correctTiles(desc,dims)
% flip 
if nargin<2
    dims = [1024 1536 251];
end
if isempty(desc)
    return;
end
desc(:,1:2) = dims(1:2)+1 - desc(:,1:2);
end