function volumeViewer_less(data, image_size, pos000)
% VOLUMEVIEWER_LESS crops 3d data array and visualize it with MATLAB built
% in funciton volumeViewer. 
% data: 3d array
% image_size: scalar or 3x1 array.
% pos000: position of the upper left front point of the image
%
%
%
if nargin < 3
    pos000 = [1,1,1];
end
if length(image_size) == 1
    image_size = round(image_size * ones([1,3]));
volumeViewer(data(pos000:pos000(1) + image_size(1), pos000(2):pos000(2) + image_size(2),...
    pos000(3):pos000(3) + image_size(3)))
end
