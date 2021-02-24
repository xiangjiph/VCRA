function [locs,xshift2D,yshift2D] = fcshift(model,order,xy,dims,locs)

if isempty(xy)
    xlocs = 1:dims(1); % dims is in (x, y, z) 
    ylocs = 1:dims(2);
    [xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
    xy = [xy1(:),xy2(:)];
end
cent_yx = squeeze(mean(model(1:2,1),3)); % p1(imaging center) in x, y direction 
p2_list = (squeeze(mean(model(1:2,2),3))); % p2(curvature) in x, y direction 
disp_xy = squeeze(mean(model(1:2,3),3)); % p3(average displacement) in x, y direction 
beta = p2_list./disp_xy.^order; % not sure the meaning of beta here
[xshift2D,yshift2D] = util.shiftxy(xy,cent_yx,beta,order,dims);
% % Visualization
% figure;
% subplot(1,2,1);
% imagesc(xshift2D);
% colorbar;
% title('Deformation field in X direction');
% subplot(1,2,2);
% imagesc(yshift2D);
% colorbar
% title('Deformation field in Y direction');
idxctrl = sub2ind(dims([2 1]),round(locs(:,2)),round(locs(:,1)));
xshift = xshift2D(idxctrl);
yshift = yshift2D(idxctrl);
locs(:,1) = locs(:,1) + xshift;
locs(:,2) = locs(:,2) + yshift;
end