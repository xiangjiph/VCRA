function plot_image(imagedata, titlename)
% Visualize matrix as image. If the matrix is float, covert to uint8 first.
% titlename: strings, can be neglected. 
% figure;
if nargin < 2
    titlename = '';
end
if isfloat(imagedata)
    imagedata = im2uint8(rescale(imagedata));
end
imshow(imagedata);
title(titlename);

end

