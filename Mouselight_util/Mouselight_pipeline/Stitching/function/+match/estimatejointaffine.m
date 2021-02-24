function [scopeparams] = estimatejointaffine(paireddescriptor,neighbors,scopeloc,params,curvemodel,old)
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

xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];

scopeparams = [];
Ntiles = size(neigs,1);
edges = cell(1,Ntiles);
for ineig = 1:Ntiles
    edge = paireddescriptor{ineig}.neigs(2:3);
    counts = paireddescriptor{ineig}.count;
    edges{ineig} = [[ineig;ineig],[edge(:)],counts(:)];
    scopeparams(ineig).imsize_um = imsize_um;
    scopeparams(ineig).dims = dims;
    
end
edges = cat(1,edges{:});
edges(any(isnan(edges),2),:)=[];
G = sparse(edges(:,1),edges(:,2),edges(:,3),Ntiles,Ntiles);
G = max(G,G');

%%
parfor_progress(Ntiles)
if old
    skipinds = any(isnan(neigs4(:,[4 5])),2);
else
    skipinds = any(isnan(neigs4(:,2:3)),2)&any(isnan(neigs4(:,4:5)),2);
end
validthis = zeros(1,Ntiles);

%%
latticeZRange = unique(scopeloc.gridix(:,3));
theselayers=latticeZRange(1:end)';
Aest_big = cell(1,max(theselayers));

parfor t = theselayers
    %%
    disp(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
    ix = (scopeloc.gridix(:,3)'==t);
    idxinlayer = find(ix);
    if (sum(ix)<1)
        disp(['No tiles found in layer #' num2str(t) '!!']);
        continue;
    end
    % apply field curvature correction to each tile in layer_t and get the
    % data matrix
    %%
    [X,sdisp] = util.createDataMatrix(params,curvemodel,paireddescriptor,scopeloc,idxinlayer);
    % AX~b: A = b\X
    Aest = sdisp/X; 
    
    res=Aest*X-sdisp;
    res = sqrt(sum(res.^2));
    Aest2 = sdisp(:,res<1)/X(:,res<1);
    
    Aest = reshape(Aest,3,3,[]);
    Aest2 = reshape(Aest2,3,3,[]);
    Aest=Aest2;
    
    
    
    
    
    %%
    Aest_big{t} = Aest;
    
end

%%
for  t = theselayers
    Aests = Aest_big{t};
    ix = (scopeloc.gridix(:,3)'==t);
    idxinlayer = find(ix);
    for ii=1:size(Aests,3)
        scopeparams(idxinlayer(ii)).affineglFC = 1000*Aests(:,:,ii); 
    end
    
end
%%

%
% %%
% parfor ineig = 1:Ntiles
%     %%
%     scopeparams(ineig).imsize_um = imsize_um;
%     scopeparams(ineig).dims = dims;
%     neiginds = find(G(ineig,:));
%     theseinds = setdiff(neiginds,paireddescriptor{ineig}.neigs(2:3));
%
%     if skipinds(ineig)
%
%     else
%         [allX,allY,allXFC,allYFC,sdisp] = deal([]);
%         siz = zeros(1); % right/below + left/above
%         stgdisp = zeros(3,1); % right/below + left/above
%         % right adjacency
%         if ~isnan(neigs(ineig,2))
%             siz(1) = size(paireddescriptor{ineig}.onx.X,1);
%             allX = [allX;paireddescriptor{ineig}.onx.X];
%             allY = [allY;paireddescriptor{ineig}.onx.Y];
%             stgdisp(:,1) = 1000*(scopeloc.loc(neigs(ineig,2),:)-scopeloc.loc(neigs(ineig,1),:));
%             sdisp = [sdisp,1000*stgdisp(:,1)*ones(1,siz(1))];
%             validthis(ineig) = 1;
%         end
%         if ~isnan(neigs(ineig,3))
%             siz(end+1) = size(paireddescriptor{ineig}.ony.X,1);
%             allX = [allX;paireddescriptor{ineig}.ony.X];
%             allY = [allY;paireddescriptor{ineig}.ony.Y];
%             stgdisp(:,end+1) = 1000*(scopeloc.loc(neigs(ineig,3),:)-scopeloc.loc(neigs(ineig,1),:));
%             sdisp = [sdisp,1000*stgdisp(:,end)*ones(1,siz(end))];
%             validthis(ineig) = 1;
%         end
%
%         if ~isempty(theseinds) & ~old
%             % check left/above
%             for ii=1:length(theseinds)
%                 if find(paireddescriptor{theseinds(ii)}.neigs==ineig)==2 % left
%                     % left adjacency
%                     ileft = neigs4(ineig,2);
%                     if ~isnan(ileft)
%                         siz(end+1) = size(paireddescriptor{ileft}.onx.X,1);
%                         allX = [allX;paireddescriptor{ileft}.onx.Y]; % flip X<->Y
%                         allY = [allY;paireddescriptor{ileft}.onx.X];
%                         stgdisp(:,end+1) = 1000*(scopeloc.loc(ileft,:)-scopeloc.loc(neigs(ineig,1),:));
%                         sdisp = [sdisp,1000*stgdisp(:,end)*ones(1,siz(end))];
%                     end
%                 else find(paireddescriptor{theseinds(ii)}.neigs==ineig)==3 % above
%                     iabove = neigs4(ineig,3);
%                     siz(end+1) = size(paireddescriptor{iabove}.ony.X,1);
%                     allX = [allX;paireddescriptor{iabove}.ony.Y]; % flip X<->Y
%                     allY = [allY;paireddescriptor{iabove}.ony.X];
%                     stgdisp(:,end+1) = 1000*(scopeloc.loc(iabove,:)-scopeloc.loc(neigs(ineig,1),:));
%                     sdisp = [sdisp,1000*stgdisp(:,end)*ones(1,siz(end))];
%                 end
%             end
%         end
%
%         if isfield(params,'beadparams') & ~isempty(params.beadparams)
%             % append bead params if exists
%             beadparamsZmatch_X = params.beadparams.allX{3};
%             beadparamsZmatch_Y = params.beadparams.allY{3};
%             beadparamsZmatchdispstage = params.beadparams.dispstage{3};
%         else
%             num = max(1,round(size(allX,1)/2));
%             beadparamsZmatch_X=ones(num,1)*[0 0 250];
%             beadparamsZmatch_Y=ones(num,1)*[0 0 0];
%             beadparamsZmatchdispstage = [ones(num,1)*[0 0 250]*1e3]';
%         end
%         siz(end+1) = size(beadparamsZmatch_X,1);
%         allX = [allX;beadparamsZmatch_X];
%         allY = [allY;beadparamsZmatch_Y];
%         sdisp = [sdisp,beadparamsZmatchdispstage];
%         allXFC = allX;
%         allYFC = allY;
%         suballX = allX+1;
%         suballY = allY+1;
%
%         % fix FC
%         [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,ineig),order,xy,dims,suballX);
%         allXFC = locs-1;
%         [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,ineig),order,xy,dims,suballY);
%         allYFC = locs-1;
%
%         Dall = (allY-allX)';
%         DallFC = (allYFC-allXFC)';
%         inds=zeros(1,sum(siz(:)));
%         siz12= [sum(siz(:,1:2),2) siz(:,3)]';
%         idxsub = [0;cumsum(siz12(:))]';
%         for ii=1:2:length(idxsub)-1
%             inds([idxsub(ii)+1:idxsub(ii+1)]) = 1;
%         end
%         dataArray.fcX = allXFC;
%         dataArray.fcY{ineig} = allXFC;
%
%         glS_=sdisp/Dall;
%         glSFC_=sdisp/DallFC;
%         scopeparams(ineig).affinegl = glS_;
%         scopeparams(ineig).affineglFC = glSFC_;
%     end
%     parfor_progress;
% end
% parfor_progress(0);
end
