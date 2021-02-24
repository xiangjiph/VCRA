function [ scopeparams_ ,paireddescriptor,curvemodel] = affineOutlierElimination( scopeloc,scopeparams,paireddescriptor,curvemodel )
%AFFINEOUTLIERELIMINATION Summary of this function goes here
%   Detailed explanation goes here

%% mean transformation
afsum = zeros(3);
Nneig = length(scopeparams);
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
%%
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

