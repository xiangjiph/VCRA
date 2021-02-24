function [paireddescriptor,curvemodel,unreliable] = fun_curvature_outlier_elimination(paireddescriptor,curvemodel,scopeloc)
%%
Nneig = size(curvemodel,3);
if 1
    %% outlier rejection based on curve fit
    unreliable = zeros(Nneig,2);
    for ineig = 1:Nneig
        if isempty(paireddescriptor{ineig}.onx.valid)
            unreliable(ineig,1) = 1;
        elseif paireddescriptor{ineig}.onx.valid==0
            unreliable(ineig,1) = 1;
        elseif isempty(paireddescriptor{ineig}.ony.valid)
            unreliable(ineig,2) = 1;
        elseif paireddescriptor{ineig}.ony.valid==0
            unreliable(ineig,2) = 1;
        end
    end
    inliers = find(~any(unreliable,2));
else
    %% outlier rejection based on median model
    validtiles = squeeze(all(curvemodel(1:2,1,:)|curvemodel(1:2,3,:)));
    [thrs,medcurvemodel] = estimateFCthreshols(curvemodel(:,:,validtiles),.99);
    medcurvemodelcents = medcurvemodel(1:2,[1 3]);
    thrcurvemodelcents = thrs(1:2,[1 3]);
    unreliable = zeros(Nneig,1);
    for ineig = 1:Nneig
        fccent = curvemodel(1:2,[1 3],ineig);
        if any(abs(fccent(:)-medcurvemodelcents(:))>thrcurvemodelcents(:))
            unreliable(ineig) = 1;
        end
    end
    inliers = find(~unreliable);
end
%%
% for every tiles estimate an affine
anchors = scopeloc.gridix(inliers,1:3);
queries = scopeloc.gridix(:,1:3);
IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]

% fill missing 
for ineig = 1:Nneig
    ianch = inliers(IDX(ineig));
    paireddescriptor{ineig}.onx = paireddescriptor{ianch}.onx;
    paireddescriptor{ineig}.ony = paireddescriptor{ianch}.ony;
    paireddescriptor{ineig}.count = [paireddescriptor{ineig}.onx.count, paireddescriptor{ineig}.ony.count];
    curvemodel(:,:,ineig) = curvemodel(:,:,ianch);
end