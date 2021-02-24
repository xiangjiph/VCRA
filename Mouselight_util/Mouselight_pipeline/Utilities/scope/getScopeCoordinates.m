function [scope] = getScopeCoordinates(inputfolder,newdash)
%GETSCOPECOORDINATES Summary of this function goes here
%
% [OUTPUTARGS] = GETSCOPECOORDINATES(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/09/21 16:33:49 $	$Revision: 0.1 $
% Copyright: HHMI 2016

args.level = 3;
args.ext = 'acquisition';
opt.seqtemp = fullfile(inputfolder,'scopeacquisitionlist.txt');
opt.inputfolder = inputfolder;
if exist(opt.seqtemp, 'file') == 2
    % load file directly
else
    args.fid = fopen(opt.seqtemp,'w');
    recdir(opt.inputfolder,args)
end

fid=fopen(opt.seqtemp,'r');
inputfiles = textscan(fid,'%s');
inputfiles = inputfiles{1};
fclose(fid);

[gridix,loc] = deal(cell(1,size(inputfiles,1)));
if newdash
    parfor_progress(size(inputfiles,1));
    parfor ifile = 1:size(inputfiles,1)
        parfor_progress;
        scvals = scopeparser(inputfiles{ifile});
        gridix{ifile} = [scvals.x scvals.y scvals.z scvals.cut_count];
        loc{ifile} = [scvals.x_mm scvals.y_mm scvals.z_mm];
    end
    parfor_progress(0);
    grids = cat(1,gridix{:});
    locs = cat(1,loc{:});
else
    parfor_progress(size(inputfiles,1));
    parfor ifile = 1:size(inputfiles,1)
        parfor_progress;
        scvals = scopeparser(inputfiles{ifile});
        loc{ifile} = [scvals.x_mm scvals.y_mm scvals.z_mm];
    end
    parfor_progress(0);
    locs = cat(1,loc{:});
    locs_ = cat(1,loc{:});
    mins = min(locs_);
    locs_ = locs_-ones(size(locs_,1),1)*mins;
    
    
    difss = diff(locs_);
    r(1) = median(difss(:,1));
    r(2) = median(difss(difss(:,2)>eps,2));
    clear x y
    [x,y] = deal(zeros(size(locs_,1),1));
    for ii=1:size(locs_)
        x(ii) = fix((locs_(ii,1)+r(1)/2)./r(1))+1;
        y(ii) = fix((locs_(ii,2)+r(2)/2)./r(2))+1;
    end
    ran = max([y x],[],1);
    oldind = 0;
    currz = 1;
    grids = zeros(size(locs_));
    for ii=1:size(locs_)
        currind = sub2ind(ran([2 1]),x(ii),y(ii));
        if currind<oldind
            % append
            currz = currz+1;
        end
        grids(ii,:) = [x(ii),y(ii),currz];
        oldind = currind;
    end
    
end

relativepaths = cell(length(inputfiles),1);
for ii=1:length(inputfiles)
    relativepaths{ii} = fileparts(inputfiles{ii}(length(inputfolder)+1:end));
end

scope.gridix = grids;
scope.loc = locs;
scope.filepath = inputfiles;
scope.relativepaths = relativepaths;
end
