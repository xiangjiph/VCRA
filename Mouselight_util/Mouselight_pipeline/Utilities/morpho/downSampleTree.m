function [outtree,selectedIdx] = downSampleTree(intree,opt)
%DOWNSAMPLETREE Summary of this function goes here
%
% [OUTPUTARGS] = DOWNSAMPLETREE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/03/29 17:01:13 $	$Revision: 0.1 $
% Copyright: HHMI 2016
[L,list] = getBranches(intree.dA);
XYZ = [intree.X intree.Y intree.Z];
R = intree.R;
D = intree.D;
%%
% down sample branches
selectedIdx = [];
for ii=1:length(L)
    %%
    set_ii = L(ii).set;
    % set_ii: +----x, includes child/excludes parent: indicies are
    % towards root
    % spacing is a function of branch length
    
    if isempty(set_ii)
        continue
    end
    
    % flip indicies so that lower indicies are close to branching
    set_ii = set_ii(end:-1:1);
    % get length of branch
    xyz = XYZ(set_ii,:);
    xyz = xyz.*(ones(size(xyz,1),1)*opt.params.voxres);
    switch opt.sampling
        case 'uni'
            p2 = sum((xyz(1:end-1,:)-xyz(2:end,:)).^2,2);
        case 'curv'
            p2 = (xyz(1:end-1,:)-xyz(2:end,:)).^2*[1;1;50];
        otherwise
    end
    dists = sqrt(p2); % resolution is 0.33 um per pixel => multiply with 3 to make it unit um
    cdist = cumsum(dists);
    
    if isempty(cdist) | cdist(end)<=1.5*opt.lengthThr %(um) : branch is a connactor and small
        % get the last index as an indicator
        indicies = set_ii(end); % keep the branch point (set flows from child towards parent)
    elseif cdist(end)>1.5*opt.lengthThr & cdist(end)<opt.largesampling
        %sample every 15(um): sampling interval for short branches
        sampling = min(1.5*opt.lengthThr,15);
        [indicies] = sampleXum(cdist,sampling,set_ii);
    elseif cdist(end)>opt.largesampling %(um)
        %sample every 50(um): sampling interval for long branches
        if 1.5*opt.lengthThr>opt.largesampling
            warning('pruning length is too big: %d compared to largesamling: %d',opt.lengthThr,opt.largesampling)
        end
        sampling = min(max(1.5*opt.lengthThr,50),opt.largesampling);
        [indicies] = sampleXum(cdist,sampling,set_ii);
    end
    selectedIdx{ii} = indicies;
end
%%
clear edges
for ii=1:length(L)
    inds = [L(ii).parentnode;selectedIdx{ii}(:)];
    to = inds(1:end-1);
    from = inds(2:end);
    edges{ii} = [from(:) to(:)]';
end
E = [edges{:}]';
%%
[idx,ia,ic] = unique(E(:));
E_ = reshape(ic,[],2);
A_ = sparse(E_(:,1),E_(:,2),1,length(idx),length(idx));
%%
outtree.dA = A_;
outtree.X = XYZ(idx,1);
outtree.Y = XYZ(idx,2);
outtree.Z = XYZ(idx,3);
outtree.R = R(idx);
outtree.D = D(idx);

end
