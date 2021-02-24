function [regpts, featmap] = loadMatchedFeatures(scopeloc,descriptorfolder,directions,checkversion,featmap)
% loadMatchedFeatures load matched features and downsample the edge voxels
% for control point generation. The downsampling is handled in this
% function, but no outlier rejection is applied 
%% Paramters
numTiles = size(scopeloc.loc,1);
downsample_opt = struct;
downsample_opt.sample_block_scale = 10;
downsample_opt.sample_block_size = [3, 3, 1] .* downsample_opt.sample_block_scale;
downsample_opt.num_sample_per_block = ceil(prod(downsample_opt.sample_block_size) * 1e-3);
downsample_opt.downsample_edgeQ = true;
downsample_opt.downsample_skelQ = false;

downsample_Q = downsample_opt.downsample_edgeQ || downsample_opt.downsample_skelQ;
%%
if nargin<3
    directions='Z'; % only check z
    checkversion = 1;
    for idx=1:numTiles
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
elseif nargin<4
    checkversion = 1;
    for idx=1:numTiles
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
elseif nargin<5
    for idx=1:numTiles
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
end
%%
tmp=cell(1,numTiles);
parfor idx = 1:numTiles
    
    matchfile = fullfile(descriptorfolder,scopeloc.relativepaths{idx},sprintf('match-%s.mat',directions));
    if exist(matchfile,'file')
        % load descriptors
        tmp{idx} = load(matchfile);
    end
    if checkversion
        matchfile2 = fullfile(descriptorfolder,scopeloc.relativepaths{idx},sprintf('match-%s-%d.mat',directions,checkversion));
        if exist(matchfile2,'file')
            % load descriptors
            tmp2 = load(matchfile2);
            if size(tmp{idx}.paireddescriptor.X,1)<size(tmp2.paireddescriptor.X,1)
                % overwrite
                tmp{idx} = tmp2;
            end
        end
    end
end
%% Down sample the paireddescriptor here
if downsample_Q
    for idx = 1 : numTiles
        if isempty(tmp{idx}) || ~isfield(tmp{idx}, 'paireddescriptor')
            continue;
        end
        tmp_paired_descriptor = tmp{idx}.paireddescriptor;
        if ~isempty(tmp_paired_descriptor.X)
%             [tmp_paired_descriptor.X, sampled_ind] = fun_uniform_sample_points_in_space(tmp_paired_descriptor.X, sample_block_size, num_sample_per_block, 'random');
%             tmp_paired_descriptor.Y = tmp_paired_descriptor.Y(sampled_ind, :);
%             tmp{idx}.paireddescriptor = tmp_paired_descriptor;
              tmp{idx}.paireddescriptor = fun_sample_matched_descriptor(tmp_paired_descriptor, downsample_opt);
        end
    end
end
%%
for idx = 1:numTiles
    featmap(idx).(genvarname(directions)) = tmp{idx};
end

% legacy variable
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
regpts = cell(1,length(featmap));for ii=1:length(regpts);if ~isfield(regpts{ii},'X');regpts{ii}.X=[];regpts{ii}.Y=[];regpts{ii}.matchrate=0;end;end
for ii=1:length(regpts)
    if isempty(featmap(ii).Z);continue;end
    regpts{ii} = featmap(ii).Z.paireddescriptor;
    regpts{ii}.neigs = [ii neighbors(ii,[4 5 7])];
end