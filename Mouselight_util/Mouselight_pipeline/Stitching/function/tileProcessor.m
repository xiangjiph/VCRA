function [descriptors,paireddescriptor,curvemodel,scopeparams] = tileProcessor(scopeloc,descriptorfolder,desc_ch,params)

checkthese = [1 4 5 7]; % 0 - right - bottom - below
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
% accumulator for steps
[paireddescriptor,R,curvemodel,scopeparams]=deal([]);

%% get tile descriptors
descriptors = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch);
sprintf('Loaded descriptors')
%% descriptor match 3785:3789
[paireddescriptor{end+1},R{end+1},curvemodel{end+1}] = match.xymatch(descriptors,neighbors(:,checkthese),scopeloc,params);
sprintf('X&Y descriptor match done')

%%
[paireddescriptor{end+1},curvemodel{end+1},unreliable] = match.curvatureOutlierElimination(paireddescriptor{end},curvemodel{end},scopeloc);
sprintf('outlier elimination done')

%%
% tile base affine
if params.singleTile
    [scopeparams{1},scopeparams{2},paireddescriptor{end+1},curvemodel{end+1}] = homographyPerTile6Neighbor(...
        params,neighbors,scopeloc,paireddescriptor{end},R,curvemodel{end});
    sprintf('per-tile affine estimation')
else
    % joint affine estimation
    [scopeparams{end+1}] = match.estimatejointaffine(paireddescriptor{end},neighbors,scopeloc,params,curvemodel{end},0);
    [scopeparams{end+1}, paireddescriptor{end+1}, curvemodel{end+1}] = match.affineOutlierElimination( scopeloc, scopeparams{end}, paireddescriptor{end},curvemodel{end},unreliable );
    sprintf('joint affine estimation')
end


