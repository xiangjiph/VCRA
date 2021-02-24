function itk_file_string = fun_vis_grid(dataset_name, stack, grid_version, vis_image, downsample_rate, overwriteQ)
% fun_vis_grid visualize the grid ( can be combined grid ) by coloring the
% space inside each box and overlaying them on the downsampled images. The
% visualization involves itk-snap. 
% If the color-coded block mask does not exist, create one and save it to
% scratch. If the color-coded mask exist, use the existing one. 
% Input: 
%   dataset_name : char array, category of the dataset
%   stack : char array, name of the image stack
%   grid_version : version of the grid to be viusalize
%   vis_image : 3D numerical array. Downsampled image of the data set
%   downsample_rate : numerical scalar. Scaling ratio, size(original
%   dataset) / size(downsampled dataset)
% Output:
%   itk_file_string: char array, filepath to the grid mask 
% Implemented by Xiang Ji on 04/08/2019

% dataset_name = 'WholeBrain';
% stack = 'mouselight_1';
% grid_version = '240_cube_combined_5';
% vis_image = DataManager.load_single_tiff(vis_image_fp);
% downsample_rate = 16;
if nargin < 4
    vis_image = [];
    downsample_rate = 1;
    overwriteQ = false;
elseif nargin < 6
    overwriteQ = true;
end
DataManager = FileManager;
itk_file_string = fullfile(DataManager.fp_quick_access, sprintf('%s_%s_%s_visualize_grid',...
    dataset_name, stack, grid_version));
itk_image_fp = sprintf('%s_image.nii', itk_file_string);
itk_mask_fp = sprintf('%s_mask.nii.gz', itk_file_string);
itk_path = DataManager.fp_itksnap_exe();
if isfile(itk_image_fp) && isfile(itk_mask_fp) && ~overwriteQ
    str_cmd = sprintf('%s -g %s -s %s &', itk_path, itk_image_fp, itk_mask_fp);
    system(str_cmd);
else
    if ~isempty(vis_image)
        image_size = size(vis_image);
        grid_c = DataManager.load_grid(dataset_name, stack, grid_version);        
        if isfield(grid_c, 'bbox_xyz_mmxx_pixel')
            all_bbox_xyz_mmxx_pixel = cat(1, grid_c.bbox_xyz_mmxx_pixel{:});
        else
            all_bbox_xyz_mmxx_pixel = cat(1, grid_c.bbox_xyz_mmxx{:});
        end
        all_bbox_xyz_mmxx_pixel_d16 = ceil(all_bbox_xyz_mmxx_pixel./downsample_rate);
        all_bbox_xyz_mmxx_pixel_d16(:, 4:6) = bsxfun(@min, all_bbox_xyz_mmxx_pixel_d16(:, 4:6), ...
            image_size);
        num_valid_bbox = size(all_bbox_xyz_mmxx_pixel_d16, 1);
        vis_mask = zeros(image_size, 'uint16');
        for idx = 1 : num_valid_bbox
            min_1 = all_bbox_xyz_mmxx_pixel_d16(idx, 1);
            min_2 = all_bbox_xyz_mmxx_pixel_d16(idx, 2);
            min_3 = all_bbox_xyz_mmxx_pixel_d16(idx, 3);
            max_1 = all_bbox_xyz_mmxx_pixel_d16(idx, 4);
            max_2 = all_bbox_xyz_mmxx_pixel_d16(idx, 5);
            max_3 = all_bbox_xyz_mmxx_pixel_d16(idx, 6);
            vis_mask(min_1:max_1, min_2:max_2, min_3:max_3) = idx;
        end
        DataManager.visualize_itksnap(vis_image, vis_mask, itk_file_string, true);
    else
        error('Missing inputs');
    end
    
end
end