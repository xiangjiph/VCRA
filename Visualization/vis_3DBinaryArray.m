function [] = vis_3DBinaryArray(binaryArray3D)
% vis_3DBinaryArray visualizes the input 3D binary array. The array is
% converted to uint8 first, smoothed with box of size 3 and visualized with
% isosurface.
% Author: Xiang Ji 10/18/2017
figure;
binaryArray3D = smooth3(uint8(binaryArray3D)*255, 'box',3);
isosurface(binaryArray3D);
daspect([1 1 1])
end
