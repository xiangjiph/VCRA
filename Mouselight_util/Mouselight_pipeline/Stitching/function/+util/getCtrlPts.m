function [toppts,xg,yg] = getCtrlPts(xpix, ypix, params,st,ed)
% calculate the locations, in pixels, of the control points (same for every
% tile).  This depends on the tile dimensions and the number of subvolumes
% to create
N = params.Ndivs;
if (xpix/N)~=round(xpix/N) || (ypix/N)~=round(ypix/N)
    error(['Image dimensions not divisible by ' num2str(N)]);
end

if nargin<4
    x_step = xpix/N;
    y_step = ypix/N;
    xg = 0:x_step-1:xpix;
    yg = 0:y_step-1:ypix;
else
    sdif = ed - st;
    x_step = sdif(1)/(N);
    y_step = sdif(2)/(N);
    xg = st(1):x_step:ed(1);
    yg = st(2):y_step:ed(2);
end
[xgrid, ygrid] = meshgrid(xg, yg);
toppts(:,1) = xgrid(:);
toppts(:,2) = ygrid(:);
toppts(:,3) = 0;
[~, IX] = sort(toppts(:,2), 'ascend');
toppts = toppts(IX,:);
end