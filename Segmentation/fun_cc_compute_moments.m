function [im_200, im_110, im_101, im_020, im_011, im_002] = fun_cc_compute_moments(cc)
% Compute moment for each connected components in three dimension. 
% Input: 
%   cc: connected component structure, output by bwconcomp
% Output:
%   im_200:
%   im_020:
%   im_002:
%   im_110:
%   im_011:
%   im_101:

num_cc = cc.NumObjects;
im_200 = zeros(num_cc,1);
im_020 = zeros(num_cc,1);
im_002 = zeros(num_cc,1);
im_110 = zeros(num_cc,1);
im_101 = zeros(num_cc,1);
im_011 = zeros(num_cc,1);
for cc_idx = 1 : num_cc
    [pixel_pos_1, pixel_pos_2, pixel_pos_3] = ind2sub(cc.ImageSize, cc.PixelIdxList{cc_idx});
    central_pos1 = mean(pixel_pos_1);
    central_pos2 = mean(pixel_pos_2);
    central_pos3 = mean(pixel_pos_3);
    pixel_pos1 = pixel_pos_1 - central_pos1;
    pixel_pos2 = pixel_pos_2 - central_pos2;
    pixel_pos3 = pixel_pos_3 - central_pos3;
    
    im_110(cc_idx) = mean(pixel_pos1 .* pixel_pos2);
    im_101(cc_idx) = mean(pixel_pos1 .* pixel_pos3);
    im_011(cc_idx) = mean(pixel_pos2 .* pixel_pos3);
    im_200(cc_idx) = mean(pixel_pos1.^2);
    im_020(cc_idx) = mean(pixel_pos2.^2);
    im_002(cc_idx) = mean(pixel_pos3.^2);
end

end