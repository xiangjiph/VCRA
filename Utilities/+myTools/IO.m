classdef IO < myTools.SystemSetUp
    % Constructor
    methods
        function obj = IO()
            obj@myTools.SystemSetUp();
        end        
    end

    methods
%% Subfolder path
        function fp = fp_constructor(obj, dataset_name, stack, datatype, file_name, makedir)
            if nargin < 6
                makedir = true;
            end
            folder_path = fullfile(obj.ROOTPATH, dataset_name, stack, datatype);
            fp = fullfile(folder_path, file_name);
            if ~isfolder(folder_path) && makedir
                    mkdir(folder_path);
            end
        end
        
        function fp = fp_Dataset(obj, dataset_name)
            fp = fullfile(obj.ROOTPATH(), dataset_name);
        end        
        
        function fp = fp_WholeBrain(obj)
            fp = obj.fp_Dataset('WholeBrain');
        end
        function fp = fp_Sample(obj)
            fp = obj.fp_Dataset('Sample');
        end
        
        function fp = fp_processed_data(obj, datasetname, stack)
            switch datasetname
                case 'WholeBrain'
                    fp = fullfile(obj.fp_Dataset(datasetname), stack, 'processed_data');
                case 'Hippocampus'
                    fp = fullfile(obj.fp_Dataset(datasetname), stack, 'processed_data');
            end
        end
        
        function fp = fp_temporary_folder(obj)
            fp = fullfile(obj.Scratch_Folder_Path, 'temp');
        end
        
        function fp = fp_quick_access(obj)
            fp = fullfile(obj.Scratch_Folder_Path, 'quick_access');
        end
        
        
%%  Raw data 
        function output = fp_raw_data_folder(obj, datasetname, stack)
            output = fullfile(obj.fp_Dataset(datasetname), stack,'raw_data');
        end
        
        function output = fp_raw_data_image(obj, datasetname,  stack, section, ext)
            if nargin<=4
                ext = '.tiff';
            end
            
            switch stack
                case 'stack1'
                    image_name = sprintf('recon_%05d%s', section, ext);
            end
               output = fullfile(obj.fp_raw_data_folder(datasetname, stack), image_name);
        end
        
        
        function output = load_single_raw_data_image(obj, datasetname, stack, section, ext)
            
            if nargin <= 4
                ext = '.tiff';
            end
            try
                fp = obj.fp_raw_data_image(datasetname, stack, section, ext);
                output = obj.load_single_tiff(fp);
            catch ME
                output = [];
                error(ME.message);
                disp('Return empty array');
            end            
        end
        
        
        function output = load_raw_data_images(obj, datasetname, stack, section_list, ext)
            if nargin <=4 
                ext = '.tiff';
            end
            fp = obj.fp_raw_data_image(datasetname, stack, section_list(1), ext);
            n_specified_slice = length(section_list);
            image_info = imfinfo(fp);
            image_size_1 = image_info(1).Height;
            image_size_2 = image_info(1).Width;
            image_bitdepth = image_info(1).BitDepth;
            switch image_bitdepth
                case 32
                    image_type = 'single';
                case 64
                    image_type = 'double';
                case 16
                    image_type = 'uint16';
                case 8
                    image_type = 'uint8';
                case 1
                    image_type = 'logical';
                otherwise
                    error('Unrecongnized image bit depth %d', image_bitdepth);
            end
            output = zeros(image_size_1, image_size_2, n_specified_slice, image_type);
            for iSection = 1 : length(section_list)
                output(:,:,iSection) = obj.load_single_raw_data_image(datasetname, stack, section_list(iSection), ext);                     
            end
        end
        
        
        function output = load_data(obj, fp)
            [~, ~, ext] = fileparts(fp);
            switch ext
                case '.mat'
                    output = load(fp);
                case {'.tiff', '.tif'}
                    output = obj.load_single_tiff(fp);
            end
        end
%%  Processed data
%%     Mask
        function output = fp_brain_mask_folder(obj,datasetname,  stack, ver)
            % stack: strings; 
            % ver: strings;
            output = fullfile(obj.fp_processed_data(datasetname, stack),'mask',ver);
        end
        
        function output = fp_brain_mask_image(obj,datasetname,  stack, ver)
            % Create filename for brain mask image tiff stack
            fn = sprintf('%s_%s_%s_mask.tiff', datasetname, stack,ver);            
            output = fullfile(obj.fp_brain_mask_folder(datasetname, stack, ver), fn);
        end
        
        function write_brain_mask(obj, data, datasetname, stack, ver)
            file_path = obj.fp_brain_mask_image(datasetname, stack, ver);
            obj.write_tiff_stack(data, file_path, 'overwrite');
        end
        
        
        function output = load_brain_mask(obj, datasetname, stack, ver, section_list)
            fp = obj.fp_brain_mask_image(datasetname, stack, ver);
            if nargin <=4
                output = obj.load_single_tiff(fp);
            else
                output = obj.load_single_tiff(fp,section_list);
            end
        end
        
        function output = fp_brain_mask_info(obj, datasetname, stack, ver)
            folder_path = obj.fp_brain_mask_folder(datasetname, stack, ver);
            fn = sprintf('%s_%s_%s_mask_info.mat', datasetname, stack,ver);
            output = fullfile(folder_path, fn);
        end
        
        
        function output = load_brain_mask_info(obj, datasetname, stack, ver)
            output = load(obj.fp_brain_mask_info(datasetname, stack, ver));
            if isfield(output, 'data')
                output = output.data;
            elseif isfield(output, 'save_data')
                output = output.save_data;
            elseif isfield(output, 'mask_info')
                output = output.mask_info;
            end
        end
        
        function write_brain_mask_info(obj, mask_info, datasetname, stack, ver) %#ok<INUSL>
            fp = obj.fp_brain_mask_info(datasetname, stack, ver);
            [foldername, ~, ~] = fileparts(fp);
            if ~isfolder(foldername)
                mkdir(foldername);
            end
            save(fp, 'mask_info');
        end       
%%     Enhanced image    
        function output = fp_enhanced_image_folder(obj, datasetname, stack)
            output = fullfile(obj.fp_processed_data(datasetname, stack), 'enhanced_images');
        end
    
        function output = fp_enhanced_image(obj, datasetname, stack, section, ext)
            % Section start from 0
            if nargin <=4 
                ext = '.tiff';
            end
            fn = sprintf('%s_%s_enhanced_%05d%s',datasetname, stack, section, ext);
            output = fullfile(obj.fp_enhanced_image_folder(datasetname, stack), fn);
        end
        
        function output = load_enhanced_image(obj, datasetname, stack, section, ext)
            if nargin <=4
                ext = '.tiff';
            end
            output = obj.load_data(obj.fp_enhanced_image(datasetname, stack, section, ext));
        end
        
        function output = fp_enhanced_image_stack_downsampled(obj, datasetname, stack, downsample_rate,info_string, ext)
            if nargin <= 4
                info_string = '';
                ext = '.tiff';
            elseif nargin <= 5
                ext = '.tiff';
            end
            fn = sprintf('%s_%s_enhanced_stack_downsampled_%dx%s%s', datasetname, stack, downsample_rate, info_string,ext);
            output = fullfile(obj.fp_enhanced_image_folder(datasetname, stack), fn);
        end
        
%%      Downsampled image        
        function output = fp_downsampled_stack_folder(obj, datasetname, stack, mask_version)
            output = fullfile(obj.fp_processed_data(datasetname, stack),'downsampled_stack',mask_version);
        end
        
        function output = fp_downsampled_stack_info_file(obj, datasetname, stack, mask_version, downsample_ratio)
            % Downsample_ratio is some number larger than 1
            file_name = sprintf('%s_%s_%s_downsampled_%dx_stack_info.mat',datasetname, stack, mask_version, downsample_ratio);
            output = fullfile(obj.fp_downsampled_stack_folder(datasetname, stack, mask_version), file_name);
        end
        
        function output = fp_downsampled_stack_file(obj, datasetname, stack, mask_version, downsample_ratio, ext, quick_accessQ)
            % Downsample_ratio is some number larger than 1
            if nargin < 6
                ext = '.tiff';
                quick_accessQ = false;
            elseif nargin < 7
                quick_accessQ = false;
            end
            file_name = sprintf('%s_%s_%s_downsampled_%dx_stack_image%s',datasetname, stack, mask_version, downsample_ratio, ext);
            if quick_accessQ
                output = fullfile(obj.fp_quick_access(), file_name);
            else
            output = fullfile(obj.fp_downsampled_stack_folder(datasetname, stack, mask_version), file_name);
            end
            
        end
        
        
        function write_downsampled_stack_info(obj, downsampled_stack_info, datasetname, stack, mask_version, downsample_ratio)
            if nargin < 3
                datasetname = downsampled_stack_info.dataset_name;
                stack = downsampled_stack_info.stack;
                mask_version = downsampled_stack_info.mask_version;
                downsample_ratio = downsampled_stack_info.downsample_ratio;
            end
            
            filepath = obj.fp_downsampled_stack_info_file(datasetname, stack, mask_version, downsample_ratio);
            [folder_path, ~,~] = fileparts(filepath);
            if ~isfolder(folder_path)
                warning('Target folder did not exist. Create new folder...\n');
                mkdir(folder_path);
            end
            save(filepath, 'downsampled_stack_info');            
        end
        
        function output = load_downsampled_stack_info(obj, datasetname, stack, mask_version, downsample_ratio)
            output = load(obj.fp_downsampled_stack_info_file(datasetname, stack, mask_version, downsample_ratio));
            output = output.downsampled_stack_info;
        end
        
        
        
        
        
%%     Grid
        function output = fp_grid_file(obj, datasetname, stack, version)
            file_name = sprintf('%s_%s_grid_%s.mat', datasetname, stack, version);
            output = fullfile(obj.fp_processed_data(datasetname, stack),'grid', file_name);
        end
        
        function output = load_grid(obj, datasetname, stack, version)
            output = load(obj.fp_grid_file(datasetname, stack, version),'grid_info');
            output = output.grid_info;
        end
        
        function write_grid_info(obj, grid_info, datasetname, stack, version) %#ok<INUSL>
            filepath = obj.fp_grid_file(datasetname, stack, version);
            save(filepath,'grid_info');
        end
%%     Block data
        function output = fp_block_data_folder(obj, datasetname, stack, version)
            output = fullfile(obj.fp_processed_data(datasetname, stack),'block_data', version);            
        end
        
        function output = fp_block_data_file(obj, datasetname, stack, version,layer, idx1, idx2)
            filename = sprintf('%s_%s_%s_block_data_%d_%d_%d.tiff', ...
                datasetname, stack, version, layer, idx1, idx2);
            output = fullfile(obj.fp_block_data_folder(datasetname, stack, version), filename);
        end
        
        function write_block_data_file(obj, data, datasetname, stack, version, layer, idx1, idx2)
            filepath = obj.fp_block_data_file(datasetname, stack, version, layer, idx1, idx2);
            [folder, ~, ~] = fileparts(filepath);
            if ~isfolder(folder)
                mkdir(folder);
            end
            obj.write_tiff_stack(data, filepath);
        end
        function output = load_block_data(obj, datasetname, stack, version, layer, idx1, idx2)
            filepath = obj.fp_block_data_file(datasetname, stack, version, layer, idx1, idx2);
            output = obj.load_single_tiff(filepath);
        end
        
        function output = fp_block_data_DICOM_foler(obj, datasetname, stack, version)
            output = fullfile(obj.fp_processed_data(datasetname, stack), 'block_data_DICOM',version);
        end
        
        function output = fp_block_data_DICOM_file(obj, datasetname, stack, version,layer, idx1, idx2)
            filename = sprintf('%s_%s_%s_block_data_%d_%d_%d.dcm', ...
                datasetname, stack, version, layer, idx1, idx2);
            output = fullfile(obj.fp_block_data_DICOM_foler(datasetname, stack, version), filename);
        end
        
        function write_block_data_DICOM_file(obj, data, datasetname, stack, version, layer, idx1, idx2)
            % write_block_data_DICOM_file save 3d array DATA to designated
            % folder. The data is saved as [idx1, idx2, color, idx_frame]
            % in DICOM file. 
            filepath = obj.fp_block_data_DICOM_file(datasetname, stack, version, layer, idx1, idx2);
            [folder, ~, ~] = fileparts(filepath);
            if ~isfolder(folder)
                warning('Target folder made.');
                mkdir(folder);
            end
            [lx, ly, lz] = size(data);
            dicomwrite(reshape(data, [lx, ly, 1, lz]), filepath);
        end
        
        function output = load_block_data_DICOM(obj, datasetname, stack, version, layer, idx1, idx2)
            filepath = obj.fp_block_data_DICOM_file(datasetname, stack, version, layer, idx1, idx2);
            output = dicomread(filepath);
        end
        
        function combined_block = load_blocks_data(obj, dataset_name, stack, grid_version, layer_list, idx_X_list, idx_Y_list)
            grid_info = obj.load_grid(dataset_name, stack, grid_version);
            num_layer = length(layer_list);
            num_X = length(idx_X_list);
            num_Y = length(idx_Y_list);

            block_cells = cell([num_X, num_Y, num_layer]);
            for layer_idx = 1 : length(layer_list)
                layer = layer_list(layer_idx);
                for x_count = 1 : length(idx_X_list)
                    idxX = idx_X_list(x_count);
                    for y_count = 1 : length(idx_Y_list)
                        idxY = idx_Y_list(y_count);
                        fprintf('Loading data from %d %d %d\n', layer, idxX, idxY);
                        idxY = idx_Y_list(y_count);
                        block_cells{x_count, y_count, layer_idx} = obj.load_block_data(dataset_name, stack, grid_version, layer, idxX, idxY);
                        block_cells{x_count, y_count, layer_idx} = block_cells{x_count, y_count, layer_idx}(1:end-grid_info.block_overlap, 1:end-grid_info.block_overlap, 1:end-grid_info.block_overlap);
                    end
                end
            end
            combined_block = cell2mat(block_cells);
        end
        
%%      ITK SNAP data
        function output = fp_itksnap_data_folder(obj, datasetname, stack, what)
            output = fullfile(obj.fp_processed_data(datasetname, stack),'itk', what);            
        end
        
        function output = fp_itksnap_file_string(obj, datasetname, stack, version, layer, idx1, idx2)
            filename = sprintf('itk_%s_%s_%s_%d_%d_%d', ...
                datasetname, stack, version, layer, idx1, idx2);
            output = fullfile(obj.fp_itksnap_data_folder(datasetname, stack, version), filename);
        end
        
        
        function output = fp_itksnap_image_file(obj, datasetname, stack, version,layer, idx1, idx2)
            output = sprintf('%s_image', obj.fp_itksnap_file_string(datasetname, stack, version, layer, idx1, idx2));
        end
        
        function output = fp_itksnap_mask_file(obj, datasetname, stack, version, layer, idx1, idx2, compressedQ)
            if nargin < 8
                compressedQ = true;
            end
            if compressedQ
                output = sprintf('%s_mask.nii.gz', obj.fp_itksnap_file_string(datasetname, stack, version, layer, idx1, idx2));
            else
                output = sprintf('%s_mask.nii', obj.fp_itksnap_file_string(datasetname, stack, version, layer, idx1, idx2)); 
            end
        end
        
        function output = load_itksnap_mask(obj, datasetname, stack, version, layer, idx1, idx2)
            output = niftiread(obj.fp_itksnap_mask_file(datasetname, stack, version, layer, idx1, idx2));
        end
        
        function write_itksnap_mask(obj, data, datasetname, stack, version, layer, idx1, idx2, add_info, compressedQ)
            itk_fp_mask = obj.fp_itksnap_file_string(datasetname, stack, version, layer, idx1, idx2);
            if nargin < 9
                add_info = '';
                compressedQ = true;
            elseif nargin < 10
                itk_fp_mask = true;
            end
            itk_fp_mask = sprintf('%s_%s', itk_fp_mask, add_info);
            
            niftiwrite(data, itk_fp_mask, 'Compressed', compressedQ);
        end
        
        function blocks_info = convert_annotation_to_block_idx(obj, dataset_name, stack, mask_version, grid_version, downsample_ratio, annotation_fn, save_infoQ)
            if nargin < 8
                save_infoQ = true;
            end
            
            if ~isfile(annotation_fn)
                annotation_fn = fullfile(obj.fp_itksnap_data_folder(dataset_name, stack, 'annotated'), ...
                    annotation_fn);
                if ~isfile(annotation_fn)
                    error('%s does not exist...', annotation_fn);
                end                
            end
            grid_info = obj.load_grid(dataset_name, stack, grid_version);
            downsample_stack_info = obj.load_downsampled_stack_info(dataset_name, stack, mask_version, downsample_ratio);
            bbox_image = niftiread(annotation_fn) > 0; % Convert uint16 to logical
            blocks_info = struct;
            blocks_info.dataset_name = dataset_name;
            blocks_info.stack = stack;
            blocks_info.annotation_path = annotation_fn;
            blocks_info.grid_info = grid_info;
            valid_bbox_count = 0;
            blocks_info.bbox_xy_idx = [];
            blocks_info.bbox_xyz_mmll = [];
            for layer_idx = 1 : grid_info.size_z
                for bbox_idx = 1 : length(grid_info.bbox_xyz_mmll{layer_idx})
                    cropped_cube = crop_bbox3(bbox_image, max(1,floor(grid_info.bbox_xyz_mmll{layer_idx}(bbox_idx,:)/downsample_stack_info.downsample_ratio)),'default');
                    %         disp(['Cube size', size(cropped_cube), '  contains: ', any(cropped_cube(:))]);
                    if any(cropped_cube(:))
                        valid_bbox_count = valid_bbox_count + 1;
                        blocks_info.bbox_xy_idx = [blocks_info.bbox_xy_idx;[layer_idx, grid_info.bbox_xy_mm_id_pos{layer_idx}(bbox_idx,:)]];
                        blocks_info.bbox_xyz_mmll = [blocks_info.bbox_xyz_mmll;grid_info.bbox_xyz_mmll{layer_idx}(bbox_idx, :)];
                    end
                end
            end
            blocks_info.valid_bbox_count = valid_bbox_count;
            if valid_bbox_count
                if save_infoQ
                    [fp, fn, ~] = fileparts(annotation_fn);
                    [~, fn, ~] = fileparts(fn);
                    save_fn = fullfile(fp, sprintf('%s_block_info.mat', fn));
                    fprintf('Block information saved to %s', save_fn);
                    save(save_fn, 'blocks_info');
                end
            else
                disp('No blocks found, please check the input information...');
            end            
        end
        
        function blocks_info = load_annotation_blocks_info(obj, dataset_name, stack, fn)
            if ~isfile(fn)
                fn = fullfile(obj.fp_itksnap_data_folder(dataset_name, stack, 'annotated'), fn);
                if ~isfile(fn)
                    error('%s does not exist...', fn);
                end                
            end
            blocks_info = load(fn);
            blocks_info = blocks_info.blocks_info;
        end
    end
%% Other Static Methods
    methods(Static)
        function write_single_tiff(data, fp, mode)
            [folder_path, ~, ~] = fileparts(fp);
            if ~isfolder(folder_path)
                warning('Folder does not exist. Create folder.');
                mkdir(folder_path)
            end
            
            if nargin < 3
                mode = 'overwrite';
            end
            imwrite(data, fp, 'WriteMode', mode, 'Compression', 'None');
        end
        
        function write_tiff_stack(data, fp, image_type,mode)
            [folder_path, ~, ~] = fileparts(fp);
            if ~isfolder(folder_path)
                warning('Folder does not exist. Create folder.');
                mkdir(folder_path)
            end
            if nargin < 3
                image_type = 'grayscale';
                mode = 'overwrite';
            elseif nargin < 4
                mode = 'overwrite';
            end
                        
            switch image_type
                case 'grayscale'
                    num_z = size(data,3);
                    imwrite(data(:,:,1), fp, 'WriteMode', mode, 'Compression', 'None');
                    if num_z > 1
                        for i_stack = 2 : num_z
                            imwrite(data(:,:,i_stack), fp, 'WriteMode', 'append', 'Compression', 'None');
                        end
                    end
                case {'color','RGB'}
                    num_z = size(data,4);
                    imwrite(data(:,:,:,1), fp, 'WriteMode', mode, 'Compression', 'None');
                    if num_z > 1
                        for i_stack = 2 : num_z
                            imwrite(data(:,:,:,i_stack), fp, 'WriteMode', 'append', 'Compression', 'None');
                        end
                    end
            end
        end
        
        
        function output = load_single_tiff(data_fp, section_list)
            % load_single_tiff loads tiff (stack) at data_fp. If
            % section_list is not given, loads all section in the stack. 
            if nargin < 2
                section_list = nan;
            end
            
            try
                n_specified_slice = length(section_list);
                image_info = imfinfo(data_fp);
                num_slice = size(image_info,1);
                image_size_1 = image_info(1).Height;
                image_size_2 = image_info(1).Width;
                image_bit_depth = image_info(1).BitDepth;
                switch image_bit_depth
                    case 32
                        image_type = 'single';
                    case 64
                        image_type = 'double';
                    case 16
                        image_type = 'uint16';
                    case 8
                        image_type = 'uint8';
                    case 1
                        image_type = 'logical'; % In matlab, logical is actually uint8
                    otherwise
                        error('Unrecongnized image bit depth %d', image_bit_depth);
                end
                
                if isnan(section_list)
                    output = zeros(image_size_1, image_size_2, num_slice, image_type);
                    for iSection  = 1 : num_slice
                        output(:,:,iSection) = imread(data_fp, iSection);
                    end
                elseif n_specified_slice == 1
                    output = imread(data_fp, section_list);
                else
                    output = zeros(image_size_1, image_size_2, n_specified_slice, image_type);
                    for iSection = 1 : n_specified_slice
                        output(:,:,iSection) = imread(data_fp, section_list(iSection));
                    end
                end
            catch ME
                error(ME.message)
            end            
        end
        
        function write_DICOM_file(fp, data)
            [folder,~,~] = fileparts(fp);
            if isfolder(folder)
                warning('Target folder made...');
                mkdir(folder);
            end
            [lx, ly, lz] = size(dadta);
            dicomwrite(reshape(data, [lx, ly, 1, lz]), fp);
        end
    end
end
