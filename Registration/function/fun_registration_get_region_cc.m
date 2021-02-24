function [subregion_cc, varargout] = fun_registration_get_region_cc(registration_str, region_id)

subregion_list_ind = fun_registration_get_subregion_list_ind(registration_str, region_id);

subregion_cc = struct;
subregion_cc.NumObjects = 2;
subregion_cc.PixelIdxList = cell(2, 1);
subregion_cc.PixelIdxList{1} = cat(2, registration_str.structure_cc_ind_hemisphere{subregion_list_ind, 1});
subregion_cc.PixelIdxList{2} = cat(2, registration_str.structure_cc_ind_hemisphere{subregion_list_ind, 2});

cc_label_mask = zeros(registration_str.array_size, 'uint8');
cc_label_mask(subregion_cc.PixelIdxList{1}) = 1;
cc_label_mask(subregion_cc.PixelIdxList{2}) = 2;

if nargout > 1
    varargout{1} = cc_label_mask;
end
end