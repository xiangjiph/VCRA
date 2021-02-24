function [xshift,yshift] = shiftxy(xy,im_cent_yx,beta,order,dims)
if nargin<4
    order = 1;
    dims = [1024 1536 251];
elseif nargin < 5
    dims = [1024 1536 251];
end
% !!!!!!!!!!!!!
% @@ TODO: is centxy always positive?? maybe need abs() here?
% !!!!!!!!!!!!!
% Why weight x and weight y are treated differently? 
if im_cent_yx(2)<=eps % no curvature
    weight_x = ((xy(:,1)-dims(1)/2).^order);
else % weightx is the distance between the point and the center
    % centxy(2) is the p1 of the fitting in +y direction, which is the position of the center in x direction. ( confusing...)
    weight_x = ((xy(:,1)-im_cent_yx(2)).^order); 
end
if im_cent_yx(1)<=eps % no curvature
    weight_y = (sign((xy(:,2) - dims(2)/2)) .* abs(xy(:,2) - dims(2)/2).^order);
else
    weight_y = (sign((xy(:,2) - im_cent_yx(1))) .* abs(xy(:,2) - im_cent_yx(1)).^order);
end
xshift = beta(1) * weight_x .* ((xy(:,2) - im_cent_yx(1)).^2);
yshift = beta(2) * weight_y .* ((xy(:,1) - im_cent_yx(2)).^2);
if nargin>4
    xshift = reshape(xshift,dims([2 1]));
    yshift = reshape(yshift,dims([2 1]));
end