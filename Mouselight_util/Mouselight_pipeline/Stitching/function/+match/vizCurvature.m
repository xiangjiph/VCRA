function vizCurvature(modd)
dims = [1024 1536 251];
xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];
[locs,xshift2D,yshift2D] = util.fcshift(modd,1,xy,dims,[1 1]);
figure, 
subplot(121)
imagesc(xshift2D)
subplot(122)
imagesc(yshift2D)