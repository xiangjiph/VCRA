function [paireddescriptor,curvemodel,scopeparams] = tileProcessor_vessel(scopeloc,descriptorfolder,desc_ch,params)

checkthese = [1 4 5 7]; % 0 - right - bottom - below
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
% accumulator for steps
[paireddescriptor,R,curvemodel,scopeparams]=deal([]);
%% descriptor match:
% Internal variable for fun_xymatch_vessel:
neigs = neighbors(:,checkthese);
% Function 
tic;
[paireddescriptor{1},R{1},curvemodel{1}] = fun_xymatch_vessel_with_IO(descriptorfolder, desc_ch, neighbors(:,checkthese), scopeloc,params);
sprintf('X&Y descriptor match done')
toc;
% save('./MouseLight_repos/Test_stitiching/output_fun_xymatch.mat', 'paireddescriptor', 'R', 'curvemodel');
% load('/home/dklab/Documents/Github/MouseLight_repos/Test_stitiching/output_fun_xymatch.mat', 'paireddescriptor', 'R', 'curvemodel');
%%
[paireddescriptor{2},curvemodel{2},unreliable] = fun_curvature_outlier_elimination_vessel(paireddescriptor{1},curvemodel{1},scopeloc);
sprintf('outlier elimination done')

%%
% tile base affine
if params.singleTile
    [scopeparams{1},scopeparams{2},paireddescriptor{3},curvemodel{3}] = fun_homographyPerTile6Neighbor_vessel(...
        params,neighbors,scopeloc,paireddescriptor{2},R,curvemodel{2});
    sprintf('per-tile affine estimation')
else
    % joint affine estimation
    [scopeparams{end+1}] = match.estimatejointaffine(paireddescriptor{end},neighbors,scopeloc,params,curvemodel{end},0);
    [scopeparams{end+1}, paireddescriptor{end+1}, curvemodel{end+1}] = match.affineOutlierElimination( scopeloc, scopeparams{end}, paireddescriptor{end},curvemodel{end},unreliable );
    sprintf('joint affine estimation')
end


