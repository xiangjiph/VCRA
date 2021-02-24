function [scopeparams,scopeparams_,paireddescriptor,curvemodel] = ...
    homographyPerTile6Neighbor(params,neighbors,scopeloc,paireddescriptor,R,curvemodel)
%HOMOGRAPHYPERTILE Summary of this function goes here
% 
% [OUTPUTARGS] = HOMOGRAPHYPERTILE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/26 10:53:18 $	$Revision: 0.1 $
% Copyright: HHMI 2017
checkthese = [1 4 5 7]; % 0 - right - bottom - below
imsize_um = params.imsize_um;
neigs = neighbors(:,checkthese);%[id -x -y +x +y -z +z] format
Nneig = size(neigs,1);
stgdisp = NaN(Nneig,3);
for ineig = 1:Nneig
    idxcent = neigs(ineig,1);
    for iadj = 1:3
        idxadj = neigs(ineig,1+iadj);
        if isnan(idxadj)
        else
            stgdisp(ineig,iadj) = 1000*(scopeloc.loc(idxadj,iadj)-scopeloc.loc(idxcent,iadj));
        end
    end
end
%%
% pix resolution based on curve model
if 0 % based on overall displacement
    xyz_umperpix_model = [stagedisplacement(:)*ones(1,size(curvemodel,3))]./squeeze(curvemodel(:,3,:));
else %based on per tile
    xyz_umperpix_model = abs(stgdisp'./squeeze(curvemodel(:,3,:)));
end
foundmatch = all(isfinite(xyz_umperpix_model(1:2,:)));
xyz_umperpix_model(1:2,foundmatch);


%% OUTLIER DETECTION
if 0
    %% based on median 
    % TODO: no need for this anymore as we do this in
    % curvatureOutlierElimination function !!
    validtiles=squeeze(all(curvemodel(1:2,1,:)|curvemodel(1:2,3,:)));
    [thrs,medcurvemodel] = estimateFCthreshols(curvemodel(:,:,validtiles),.99);
    medcurvemodelcents = medcurvemodel(1:2,[1 3]);
    thrcurvemodelcents = thrs(1:2,[1 3]);
    count = zeros(Nneig,2);
    for ineig = 1:Nneig
        [count(ineig,:)] = paireddescriptor{ineig}.count;
    end
    unreliable = zeros(Nneig,1);
    for ineig = 1:Nneig
        fccent = curvemodel(1:2,[1 3],ineig);
        if any(abs(fccent(:)-medcurvemodelcents(:))>thrcurvemodelcents(:))
            unreliable(ineig) = 1;
        end
    end
    %%
    inliers = find(~unreliable);
    % for every tiles estimate an affine
    anchors = scopeloc.gridix(inliers,1:3);
    queries = scopeloc.gridix(:,1:3);
    IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
    
    %%
    % fill missing
    for ineig = 1:Nneig
        ianch = inliers(IDX(ineig));
        paireddescriptor{ineig}.onx = paireddescriptor{ianch}.onx;
        paireddescriptor{ineig}.ony = paireddescriptor{ianch}.ony;
        paireddescriptor{ineig}.count = [size(paireddescriptor{ineig}.onx.X,1) size(paireddescriptor{ineig}.ony.X,1)];
        curvemodel(:,:,ineig) = curvemodel(:,:,ianch);
    end
else
    validtiles=squeeze(all(curvemodel(1:2,1,:)|curvemodel(1:2,3,:)));
    % rejection based on projected points
end

%%
params.order = 1;
[scopeparams,validthis] = util.estimateaffine(paireddescriptor,neighbors,scopeloc,params,curvemodel,0);

%% affine outlier: TODO, reject based on tile corner pixel shift magnitudes!!
%mean transformation
afsum = zeros(3);
iter=0;
for ineig = 1:Nneig
    af = scopeparams(ineig).affineglFC;
    if isempty(af)
    else
       afsum = afsum+af;
       iter=iter+1;
    end
end
aff = afsum/iter;

%%
scopeparams_ = scopeparams;
reliable = zeros(Nneig,1);
for ineig = 1:Nneig
    if isempty(scopeparams_(ineig).affineglFC) | norm(scopeparams_(ineig).affineglFC-aff)/norm(aff)*100>1
        %not reliable
    else
        reliable(ineig) = 1;
    end
end

inliers = find(reliable(:));
% for every tiles estimate an affine
anchors = scopeloc.gridix(inliers,1:3);
queries = scopeloc.gridix(:,1:3);
IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
% fill missing 
for ineig = 1:Nneig
    ianch = inliers(IDX(ineig));
    paireddescriptor{ineig}.onx = paireddescriptor{ianch}.onx;
    paireddescriptor{ineig}.ony = paireddescriptor{ianch}.ony;
    curvemodel(:,:,ineig) = curvemodel(:,:,ianch);
    scopeparams_(ineig) = scopeparams(ianch);
end


end
