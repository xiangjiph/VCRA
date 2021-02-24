function [regpts,featmap] = loadMatchedFeatures(scopeloc,descriptorfolder,directions,checkversion,featmap)
% get list file
numTiles = size(scopeloc.loc,1);

if nargin<3
    directions='Z'; % only check z
    checkversion = 0;
    for idx=1:numTiles
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
elseif nargin<4
    checkversion = 0;
    for idx=1:numTiles
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
elseif nargin<5
    for idx=1:numTiles
        [featmap(idx).X,featmap(idx).Y,featmap(idx).Z] = deal([]);
    end
end
% myfile = dir(fullfile(descriptorfolder,'list*files'));
% fid=fopen(fullfile(descriptorfolder,myfile.name),'r');
% inputfiles = textscan(fid,'%s');inputfiles = inputfiles{1};fclose(fid);
% folders = unique(cellfun(@fileparts,inputfiles,'UniformOutput',false),'stable');

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

