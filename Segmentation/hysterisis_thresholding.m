function mask = hysterisis_thresholding(image,th_low,th_high)
% HYSTERISIS_THRESHOLDING computes the mask by hysterisis thresholding.
% Input:
%   th1: lower threshold 
%   th2: higher threshold
%   Both threshold can be a number of an array of the same size as the
%   image.
% Output: 
%   mask: logical array 
% Note: No need to implement on GPU. For 512 ^ 3 block, it takes CPU ~1.5
% seconds while GPU needs ~10 seconds. The main reason might be memory
% access at the last step of this implementation.

if th_low(1) > th_high(1)
    warning('The input lower threshold is larger than the input higher threshold. Swapped automatically');
    tmp = th_low;
    th_low = th_high;
    th_high = tmp;
    clear tmp;
end
mask = zeros(size(image), 'logical');
mask_low = image >= th_low;
mask_high = image >= th_high;
if isa(mask_low, 'gpuArray')
    mask_low = gather(mask_low);
    mask_high = gather(mask_high);
end

low_cc = bwconncomp(mask_low, 8);
for cc_idx = 1 : low_cc.NumObjects
    if sum(mask_high(low_cc.PixelIdxList{cc_idx}))>0
        mask(low_cc.PixelIdxList{cc_idx}) = true;
    end
end


end