function [paireddescriptor,curvemodel,unreliable] = fun_curvature_outlier_elimination_vessel(paireddescriptor,curvemodel,scopeloc)
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
        end
        
        if isempty(paireddescriptor{ineig}.ony.valid)
            unreliable(ineig,2) = 1;
        elseif paireddescriptor{ineig}.ony.valid==0
            unreliable(ineig,2) = 1;
        end
    end
    unreliable_x = unreliable(:,1);
    unreliable_y = unreliable(:,2);
    inliers_x = find(~unreliable_x);
    inliers_y = find(~unreliable_y);
    inliers_xy = find(~any(unreliable,2));
else
    %% outlier rejection based on median model
    validtiles = squeeze(all(curvemodel(1:2,1,:)|curvemodel(1:2,3,:)));
    [thrs,medcurvemodel] = estimateFCthreshols(curvemodel(:,:,validtiles),.99);
    medcurvemodelcents = medcurvemodel(1:2,[1 3]);
    thrcurvemodelcents = thrs(1:2,[1 3]);
    unreliable = zeros(Nneig,1);
    for ineig = 1:Nneig
        fccent = curvemodel(1:2,[1 3],ineig);
        if any(abs(fccent(:) - medcurvemodelcents(:)) > thrcurvemodelcents(:))
            unreliable(ineig) = 1;
        end
    end
    inliers = find(~unreliable);
end
%%
% For each tile, each direction, find the nearest neighbor
% X direction 
anchors_x = scopeloc.gridix(inliers_x,1:3);
queries_x = scopeloc.gridix(:,1:3);
knn_IDX_x = knnsearch(anchors_x,queries_x,'K',1,'distance',@distfun);%W=[1 1 100000]
% Y direction 
anchors_y = scopeloc.gridix(inliers_y,1:3);
queries_y = scopeloc.gridix(:,1:3);
knn_IDX_y = knnsearch(anchors_y,queries_y,'K',1,'distance',@distfun);%W=[1 1 100000]
% for every tiles estimate an affine
anchors = scopeloc.gridix(inliers_xy,1:3);
queries = scopeloc.gridix(:,1:3);
IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]

% fill missing 
for ineig = 1 : Nneig
%     if unreliable_x(ineig)
%         ianch_x = inliers_x(knn_IDX_x(ineig));
%         paireddescriptor{ineig}.onx = paireddescriptor{ianch_x}.onx;
%         curvemodel(1,:,ineig) = curvemodel(1,:,ianch_x);
%     end
%     if unreliable_y(ineig)
%         ianch_y = inliers_y(knn_IDX_y(ineig));
%         paireddescriptor{ineig}.ony = paireddescriptor{ianch_y}.ony;
%         curvemodel(2,:,ineig) = curvemodel(2,:,ianch_y);
%     end    
    ianch = inliers_xy(IDX(ineig));
    paireddescriptor{ineig}.onx = paireddescriptor{ianch}.onx;
    paireddescriptor{ineig}.ony = paireddescriptor{ianch}.ony;
    curvemodel(:,:,ineig) = curvemodel(:,:,ianch);
    paireddescriptor{ineig}.count = [size(paireddescriptor{ineig}.onx.X,1) size(paireddescriptor{ineig}.ony.X,1)];
end
%% Check the tile that are replaced
% scopeloc_min = min(scopeloc.gridix(:,1:3), [], 1);
% inliers_grid_sub = (anchors - scopeloc_min + 1);
% inliers_grid_ind = sub2ind([7,7,7], inliers_grid_sub(:,1), inliers_grid_sub(:,2), ...
%     inliers_grid_sub(:,3));
% inliers_grid_mask = false([7,7,7]);
% inliers_grid_mask(inliers_grid_ind) = true;
% implay(inliers_grid_mask)
% test_grid_sub = [6,3,6];
% test_grid_sub_global = test_grid_sub + scopeloc_min - 1;
% test_list_idx = find( (scopeloc.gridix(:,1) == test_grid_sub_global(1)) & ...
%     (scopeloc.gridix(:,2) == test_grid_sub_global(2)) & ...
%     (scopeloc.gridix(:,3) == test_grid_sub_global(3)));
% test_paired_descriptor = paireddescriptor{test_list_idx};
%%
% % scatter3(skl_sub(:,1), skl_sub(:,2), skl_sub(:,3))
% figure;
% plot_sub_1 = test_paired_descriptor.onx.X_skl;
% plot_sub_2 = test_paired_descriptor.onx.Y_skl;
% scatter3(plot_sub_1(:,2), plot_sub_1(:,1), plot_sub_1(:,3));
% hold on 
% scatter3(plot_sub_2(:,2), plot_sub_2(:,1), plot_sub_2(:,3));
end