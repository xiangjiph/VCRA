classdef FileManager

%% Properties
    properties
        computer_name;
        ROOT_PATH;
        SCRIPT_PATH;
        EXTERNAL_LIB_ROOT_PATH;
        SCRATCH_ROOT_PATH;
        SERVER_ROOT_PATH;
        SERVER_SCRIPT_PATH
        opt_load_from_scratch_first_Q logical = false;
        
    end
%% Directory settings
   methods
       function obj = FileManager(host_name, load_from_scratch_Q, ~)
           if nargin < 1
               host_name = [];
               load_from_scratch_Q = false;
%                write_to_sratch_Q = false;
           elseif nargin < 2
               load_from_scratch_Q = false;
           end
           obj.opt_load_from_scratch_first_Q = load_from_scratch_Q;
           if isempty(host_name)
               [~, host_name] = system('hostname');
               host_name = host_name(1:end-1);
           end
           obj.computer_name = host_name;
           switch host_name
               case {'MACHINE 1'}
                   obj.ROOT_PATH = 'ABSOLUTE PATH TO THE ROOT OF THE DATA FOLDER';
                   obj.SCRIPT_PATH = 'ABSOLUTE PATH TO THE ROOT OF THE CODE';
                   obj.EXTERNAL_LIB_ROOT_PATH = 'ABSOLUTE PATH TO THE ROOT OF THE EXTERNAL LIBRARIES (e.g. ITK-SNAP)';
                   obj.SCRATCH_ROOT_PATH = 'ABSOLUTE PATH TO THE SCRATCH DISK FOLDER';
                   obj.SERVER_ROOT_PATH = 'ABSOLUTE PATH TO THE SERVER';  
                   obj.SERVER_SCRIPT_PATH = 'ABSOLUTE PATH TO THE SHARED CODE ON THE SERVER';  
               case 'MACHINE 2'
                   obj.ROOT_PATH = 'ABSOLUTE PATH TO THE ROOT OF THE DATA FOLDER';
                   obj.SCRIPT_PATH = 'ABSOLUTE PATH TO THE ROOT OF THE CODE';
                   obj.SCRATCH_ROOT_PATH = 'ABSOLUTE PATH TO THE SCRATCH DISK FOLDER';
                   obj.SERVER_ROOT_PATH = 'ABSOLUTE PATH TO THE SERVER';                   
               otherwise
                   % Use current directory
                   file_path_parts = strsplit(pwd, filesep);
                   obj.ROOT_PATH = fullfile(file_path_parts{1}, 'Vessel');
           end
           
       end
   end
%% Methods for getting file paths    
    methods 
%%  General
        function fp = fp_constructor(obj, dataset_name, stack, datatype, file_name, makedir)
            if nargin < 6
                makedir = true;
            end
            folder_path = fullfile(obj.ROOT_PATH, dataset_name, stack, datatype);
            fp = fullfile(folder_path, file_name);
            if ~isfolder(folder_path) && makedir
                    mkdir(folder_path);
            end
        end
        
        function fp = fp_Dataset(obj, dataset_name)
            fp = fullfile(obj.ROOT_PATH, dataset_name);
        end        
        
        function fp = fp_WholeBrain(obj)
            fp = obj.fp_Dataset('WholeBrain');
        end
        function fp = fp_Sample(obj)
            fp = obj.fp_Dataset('Sample');
        end
        
        function fp = fp_processed_data(obj, dataset_name, stack)
            switch dataset_name
                case 'WholeBrain'
                    fp = fullfile(obj.fp_Dataset(dataset_name), stack, 'processed_data');
                otherwise
                    fp = fullfile(obj.fp_Dataset(dataset_name), stack, 'processed_data');
            end
        end
        
        function fp = fp_temporary_folder(obj)
            fp = fullfile(obj.SCRATCH_ROOT_PATH, 'temp');
        end
        
        function fp = fp_quick_access(obj)
            fp = fullfile(obj.SCRATCH_ROOT_PATH, 'quick_access');
        end
        
        function fp = fp_file_exist_locally(obj, fp)
            if ~isfile(fp)
                switch obj.computer_name
                    case {'MACHINE1', 'OTHER NAME OF MACHINE 1'}
                        fp = strrep(fp, obj.ROOT_PATH, obj.SERVER_ROOT_PATH);
                        if ~isfile(fp)
                            warning('File exist on neither local disk nor server');
                        end
                    case 'MACHINE 2'
                        % 1. Windows linux subsystem is required to run
                        % wsl. See the following link for reference:
                        % https://docs.microsoft.com/en-us/windows/wsl/interop
                        
                        % 2. A public key must be set up for ssh before rsync.
                        % to copy files from the server. To set up the
                        % public key, see: https://www.tecmint.com/ssh-passwordless-login-using-ssh-keygen-in-5-easy-steps/
                        % Notice that when creating the ras file, do NOT
                        % rename the file, otherwise ssh will still ask for
                        % password
                        disp('File does not exist locally. Try to download from the server');
                        fp_wsl = strrep(fp, filesep, '/');
                        fp_wsl = strrep(fp_wsl, 'C:', '/mnt/c');
                        [folder_path, ~, ~] = fileparts(fp);
                        if ~isfolder(folder_path)
                            mkdir(folder_path);
                        end
                        fp_remote = strrep(fp, obj.ROOT_PATH, obj.SERVER_ROOT_PATH);
                        fp_remote = strrep(fp_remote, filesep, '/');
                        fp_remote = sprintf('MACHINE 1:%s', fp_remote);
                        command_str = sprintf('wsl rsync -arv %s %s', fp_remote, fp_wsl);
                        system(command_str);
                        if ~isfile(fp)
                            warning('File does not exist on local machien nor remote workstation');
                        end
                    otherwise
                        error('Unrecognized machine');
                end
            end
        end
        
        function fp = fp_file_exist_in_scratch(obj, fp)
           fp_scratch = strrep(fp, obj.ROOT_PATH, obj.SCRATCH_ROOT_PATH);
           if isfile(fp_scratch)
               fp = fp_scratch;
           end            
        end        
%%  Raw data 
        function output = fp_raw_data_folder(obj, dataset_name, stack)
            output = fullfile(obj.fp_Dataset(dataset_name), stack,'raw_data');
        end
%%  Processed data
%%     Mask
        function output = fp_mask_folder(obj,dataset_name,  stack, ver)
            % stack: strings; 
            % ver: strings;
            output = fullfile(obj.fp_processed_data(dataset_name, stack),'mask',ver);
        end
        
        function output = fp_brain_mask_image(obj,dataset_name,  stack, ver)
            % Create filename for brain mask image tiff stack
            fn = sprintf('%s_%s_%s_mask.tiff', dataset_name, stack,ver);            
            output = fullfile(obj.fp_mask_folder(dataset_name, stack, ver), fn);
        end
        
        function output = fp_block_mask_file(obj, dataset_name, stack, version, idx1, idx2, layer, ext)
            if nargin < 8
                ext = '.mat';
            end
            
            filename = sprintf('%s_%s_%s_block_data_%d_%d_%d%s', ...
                dataset_name, stack, version, idx1, idx2, layer, ext);
            output = fullfile(obj.fp_mask_folder(dataset_name, stack, version), filename);
        end
        
        function write_block_mask(obj, data, dataset_name, stack, version, idx1, idx2, layer)
            if isa(data, 'numeric')
                ext = '.tiff';
            elseif isa(data, 'struct')
                ext = '.mat';
            else
                error('Unrecognized data type');
            end 
            file_path = obj.fp_block_mask_file(dataset_name, stack, version, idx1, idx2, layer, ext);
            [folder_path, ~, ~] = fileparts(file_path);
            if ~isfolder(folder_path)
                warning('Target folder does not exist. Create folder.');
                mkdir(folder_path);
            end
            
            if isa(data, 'numeric')
                obj.write_tiff_stack(data, file_path, 'grayscale', 'overwrite');
            elseif isa(data, 'struct')
                save(file_path, '-struct', 'data');
            end
        end
        
        function output = load_block_mask(obj, dataset_name, stack, version, idx1, idx2, layer, ext)
            if nargin < 8
                ext = '.mat';
            end
            fp = obj.fp_block_mask_file(dataset_name, stack, version, idx1, idx2, layer, ext);
            output = obj.load_data(fp);
        end
        
        function combined_block = load_blocks_files(obj, what, dataset_name, stack, grid_version, idx_1_list, idx_2_list, layer_list, dtype, data_version)
            % One potential problem for this function is that the size of
            % the empty block might not be grid_info.block_size, if these
            % empty blocks are on the boundary of the whole data set. TO be
            % fixed later. 
            if nargin < 9
                switch what
                    case 'mask'
                        dtype = 'logical';
                    case {'image', 'data'}
                        dtype = 'uint16';
                    case {'dt', 'distance_transform', 'skl', 'skel', 'skeleton'}
                        dtype = 'single';
                end
                data_version = grid_version;
            elseif nargin < 10
                data_version = grid_version;
            end
            grid_info = obj.load_grid(dataset_name, stack, grid_version);
            
            if isfield(grid_info, 'data_type') && any(strcmp(what, {'image', 'data'}))
                dtype = grid_info.data_type;
            end 
            block_overlap = grid_info.block_overlap;
            half_block_overlap = block_overlap/2;
            num_layer = length(layer_list);
            num_X = length(idx_1_list);
            num_Y = length(idx_2_list);

            block_cells = cell([num_X, num_Y, num_layer]);
            for layer_idx = 1 : num_layer
                layer = layer_list(layer_idx);
                for x_count = 1 : num_X
                    idx1 = idx_1_list(x_count);
                    for y_count = 1 : num_Y
                        idx2 = idx_2_list(y_count);
                        if grid_info.bbox_xy_label_mat{layer}(idx1, idx2)
%                             fprintf('Loading data from %d %d %d\n', idx1, idx2, layer);
                            switch what
                                case 'mask'
                                    block_data = obj.load_block_mask(dataset_name, stack, data_version, idx1, idx2, layer);
                                    if isa(block_data, 'struct')
                                        tmp = block_data;
                                        block_data = false(tmp.block_size);
                                        block_data(tmp.ind) = true;
                                    end
                                case {'image', 'data'}
                                    block_data = obj.load_block_data(dataset_name, stack, grid_version, idx1, idx2, layer);
                                case {'dt', 'distance_transform'}
                                    dt_str = obj.load_block_dt(dataset_name, stack, grid_version, idx1, idx2, layer);
                                    if isa(dt_str, 'struct')
                                        block_data = zeros([grid_info.grid2D.ll(sub2ind(grid_info.grid2D.size, idx1,idx2),:), grid_info.gridZ.ll(layer)], dtype);
                                        block_data(dt_str.ind) = dt_str.dt;
                                    elseif isa(dt_str, 'single')
                                        block_data = dt_str;
                                    end
                                case {'skl', 'skel', 'skeleton'}
                                    skl_str = obj.load_block_skl(dataset_name, stack, data_version, idx1, idx2, layer);
                                    block_data = zeros(skl_str.block_size, 'single');
                                    block_data(skl_str.ind) = skl_str.r;
                            end
                        else
%                             warning('Block (%d, %d, %d) does not exist. Use zeros array.', idx1, idx2, layer);
                            if ~isfield(grid_info, 'grid2D')
                                block_data = zeros(grid_info.block_size, dtype);
                            else
                                block_data = zeros([grid_info.grid2D.ll(sub2ind(grid_info.grid2D.size, idx1,idx2),:), grid_info.gridZ.ll(layer)], dtype);
                            end
                        end
                        
                        if layer_idx ~= num_layer
                            block_data = block_data(:,:, 1:end-half_block_overlap);
                        end
                        if layer_idx ~= 1
                            block_data = block_data(:,:, (1 + half_block_overlap):end);
                        end
                        if x_count ~= num_X
                            block_data = block_data(1:end-half_block_overlap, :, :);
                        end
                        if x_count ~= 1
                            block_data = block_data((1 + half_block_overlap):end, :, :);
                        end
                        if y_count ~= num_Y
                            block_data = block_data(:,1:end-half_block_overlap,:);
                        end
                        if y_count ~= 1
                            block_data = block_data(:,(1 + half_block_overlap):end,:);
                        end
                        block_cells{x_count, y_count, layer_idx} = block_data;
                    end
                end
            end
            combined_block = cell2mat(block_cells);
        end
        
        function write_brain_mask(obj, data, dataset_name, stack, ver)
            file_path = obj.fp_brain_mask_image(dataset_name, stack, ver);
            obj.write_tiff_stack(data, file_path);
        end
        
        
        function output = load_brain_mask(obj, dataset_name, stack, ver, section_list)
            fp = obj.fp_brain_mask_image(dataset_name, stack, ver);
            if nargin <=4
                output = obj.load_data(fp);
            else
                fp = obj.fp_file_exist_in_scratch(fp);
                fp = obj.fp_file_exist_locally(fp);
                output = obj.load_single_tiff(fp, section_list);
            end
        end
        
        function output = fp_mask_info(obj, dataset_name, stack, ver)
            folder_path = obj.fp_mask_folder(dataset_name, stack, ver);
            fn = sprintf('%s_%s_%s_mask_info.mat', dataset_name, stack,ver);
            output = fullfile(folder_path, fn);
        end
        
        
        function output = load_mask_info(obj, dataset_name, stack, ver)
            fp = obj.fp_mask_info(dataset_name, stack, ver);
            output = obj.load_data(fp);
            if isfield(output, 'data')
                output = output.data;
            elseif isfield(output, 'save_data')
                output = output.save_data;
            elseif isfield(output, 'mask_info')
                output = output.mask_info;                
            end
        end
        
        function write_mask_info(obj, mask_info, dataset_name, stack, ver)
            fp = obj.fp_mask_info(dataset_name, stack, ver);
            [foldername, ~, ~] = fileparts(fp);
            if ~isfolder(foldername)
                mkdir(foldername);
            end
            save(fp, '-struct','mask_info');
        end      
        
        function output = fp_analysis_data_folder(obj, dataset_name, stack, subfolder_name)
            if nargin < 4
                subfolder_name = [];
            end
           output = fullfile(obj.fp_processed_data(dataset_name, stack), 'analysis_data', subfolder_name);
        end
        
        function output = fp_analysis_data_file(obj, dataset_name, stack, what, subfolder_name)
            if nargin < 5
                subfolder_name = [];
            end
            output = fullfile(obj.fp_analysis_data_folder(dataset_name, stack, subfolder_name), what);
        end
        
        function output = fp_analysis_data_file_in_grid(obj, subfolder_name, ...
                dataset_name, stack, grid_version, grid_label)            
            output = fullfile(obj.fp_analysis_data_folder(dataset_name, stack, subfolder_name), ...
                sprintf('%s_%s_%s_%d_%s.mat', dataset_name, stack, grid_version, ...
                grid_label, subfolder_name));
        end
        
        function output_fp = write_analysis_data_in_grid(obj, data, subfolder_name, dataset_name, ...
                stack, grid_version, grid_label)
            output_fp = obj.fp_analysis_data_file_in_grid(subfolder_name, ...
                dataset_name, stack, grid_version, grid_label);
            obj.write_data(output_fp, data);
        end
        
        function output = load_analysis_data_in_grid(obj, subfolder_name, dataset_name, ...
                stack, grid_version, grid_label)
            fp = obj.fp_analysis_data_file_in_grid(subfolder_name, ...
                dataset_name, stack, grid_version, grid_label);
            output = obj.load_data(fp);
        end
        
        function output = load_analysis_data(obj, dataset_name, stack, what, subfolder)
           if nargin < 5
               subfolder = [];
           end
           fp = obj.fp_analysis_data_file(dataset_name, stack, what, subfolder);
           output = obj.load_data(fp);
        end
        
        function exit_code = write_analysis_data(obj, data, dataset_name, stack, what, subfolder)
            if nargin < 6
                subfolder = [];
            end
            fp = obj.fp_analysis_data_file(dataset_name, stack, what, subfolder);
            [folder, ~, ~] = fileparts(fp);
            if ~isfolder(folder)
                mkdir(folder);
            end
            obj.write_data(fp, data);
            exit_code = 0;
        end
%%      Downsampled image        
        function output = fp_downsampled_stack_folder(obj, dataset_name, stack, mask_version)
            output = fullfile(obj.fp_processed_data(dataset_name, stack),'downsampled_stack',mask_version);
        end
        
        function output = fp_downsampled_stack_info_file(obj, dataset_name, stack, mask_version, downsample_ratio)
            % Downsample_ratio is some number larger than 1
            file_name = sprintf('%s_%s_%s_downsampled_%dx_stack_info.mat',dataset_name, stack, mask_version, downsample_ratio);
            output = fullfile(obj.fp_downsampled_stack_folder(dataset_name, stack, mask_version), file_name);
        end
        
        function output = fp_downsampled_stack_file(obj, dataset_name, stack, mask_version, downsample_ratio, ext, quick_accessQ)
            % Downsample_ratio is some number larger than 1
            if nargin < 6
                ext = '.tiff';
                quick_accessQ = false;
            elseif nargin < 7
                quick_accessQ = false;
            end
            file_name = sprintf('%s_%s_%s_downsampled_%dx_stack_image%s',dataset_name, stack, mask_version, downsample_ratio, ext);
            if quick_accessQ
                output = fullfile(obj.fp_quick_access(), file_name);
            else
            output = fullfile(obj.fp_downsampled_stack_folder(dataset_name, stack, mask_version), file_name);
            end
        end
        
        function write_downsampled_stack_info(obj, downsampled_stack_info, dataset_name, stack, mask_version, downsample_ratio)
            if nargin < 3
                dataset_name = downsampled_stack_info.dataset_name;
                stack = downsampled_stack_info.stack;
                mask_version = downsampled_stack_info.mask_version;
                downsample_ratio = downsampled_stack_info.downsample_ratio;
            end
            
            filepath = obj.fp_downsampled_stack_info_file(dataset_name, stack, mask_version, downsample_ratio);
            [folder_path, ~,~] = fileparts(filepath);
            if ~isfolder(folder_path)
                warning('Target folder did not exist. Create new folder...\n');
                mkdir(folder_path);
            end
            save(filepath,'-struct','downsampled_stack_info');            
        end
        
        function output = load_downsampled_stack_info(obj, dataset_name, stack, mask_version, downsample_ratio)
            output = obj.load_data(obj.fp_downsampled_stack_info_file(dataset_name, stack, mask_version, downsample_ratio));
            if isfield(output, 'downsampled_stack_info')
                output = output.downsampled_stack_info;
            end
        end
%%      Thumbnail
        function output = fp_thumbnail_folder(obj, dataset_name, stack)
            output = fullfile(obj.fp_processed_data(dataset_name, stack), 'thumbnail');
        end
        
        function output = fp_thumbnail_file(obj, dataset_name, stack, section, ext)
            if nargin < 5
                ext = '.tif';
            end
            file_name = sprintf('%05d%s', section, ext);
            output = fullfile(obj.fp_thumbnail_folder(dataset_name, stack), file_name);
        end
        
        function output = load_thumbnail(obj, dataset_name, stack, section)
            output = obj.load_data(obj.fp_file_exist_locally(obj.fp_thumbnail_file(dataset_name, stack, section)));            
        end
%%     Grid
        function output = fp_grid_file(obj, dataset_name, stack, version)
            file_name = sprintf('%s_%s_grid_%s.mat', dataset_name, stack, version);
            output = fullfile(obj.fp_processed_data(dataset_name, stack),'grid', file_name);
        end
        
        function output = load_grid(obj, dataset_name, stack, version)
            try 
                output = obj.load_data(obj.fp_grid_file(dataset_name, stack, version));    
            catch
                output = obj.load_data(obj.fp_grid_file(dataset_name, stack, version),'grid_info');
            end
            if isfield(output, 'grid_info')
                output = output.grid_info;
            end
            
        end
        
        function write_grid_info(obj, grid_info, dataset_name, stack, version) 
            filepath = obj.fp_grid_file(dataset_name, stack, version);
            [folder_path, ~, ~] = fileparts(filepath);
            if ~isfolder(folder_path)
                warning('Target folder does not exist. Create a new one.');
                mkdir(folder_path);
            end
            save(filepath,'-struct', 'grid_info');
        end
        
        function grid_coordinate_target = convert_between_multiresolution_grids(grid_ori, target_downsample_rate, layer_ori, ...
                idx1_ori, idx2_ori, row_ori, col_ori, sec_ori)
            pos_1x = (grid_ori.bbox_xyz_mmll{layer_ori}(grid_ori.bbox_xy_linear_idx_mat{layer_ori}(idx1_ori, idx2_ori),1:3) + [row_ori, col_ori, sec_ori])* ...
                round(grid_ori.downsample_rate/target_downsample_rate) ;
            pos_target = pos_1x ./ target_downsample_rate;
            block_size_target = grid_ori.block_size_1x - ceil(grid_ori.block_overlap_1x * (1 - 1/target_downsample_rate));
            block_overlap_target = grid_ori.block_overlap_1x/target_downsample_rate;
            block_spacing_target = block_size_target - block_overlap_target;
            idx_pos_target = ceil(pos_target ./ (block_spacing_target));
            idx_pos_target = circshift(idx_pos_target, 1);
            pos_target_in_block = mod(pos_target, block_spacing_target);
            grid_coordinate_target = [idx_pos_target, pos_target_in_block];
        end
%%      Metadata file
        function output = fp_metadata_folder(obj, dataset_name, stack)
            output = fullfile(obj.fp_processed_data(dataset_name, stack), 'metadata');
            if ~isfolder(output)
                warning('Target folder does not exist. Create one');
                mkdir(output);
            end
            
        end
        
        function output = fp_metadata_file(obj, dataset_name, stack, what, ext)
            if nargin < 5
                ext = '.mat';
            end
            [file_folder, file_name, file_ext] = fileparts(what);
            if ~isempty(file_ext)
                ext = file_ext;
                what = fullfile(file_folder, file_name);
            end
            output = fullfile(obj.fp_metadata_folder(dataset_name, stack),...
                sprintf('%s%s', what, ext));
        end
        
        function output = load_metadata(obj, dataset_name, stack, what, ext)
            if nargin < 5
                ext = '.mat';
            end
            
            fp = obj.fp_metadata_file(dataset_name, stack, what, ext);
            output = obj.load_data(fp);
        end        
        
        function output = fp_training_data(obj, dataset_name, stack, name_string)
            fn = sprintf('Training_data_%s.mat', name_string);
            output = fullfile( obj.fp_metadata_folder(dataset_name, stack), fn);
        end
        
        function output = fp_classifier(obj, dataset_name, stack, name_string)
            fn = sprintf('Classifier_%s.mat', name_string);
            output = fullfile( obj.fp_metadata_folder(dataset_name, stack), fn);
        end
                        
        function output = fp_registration_data(obj, dataset_name, stack, reg_ver)
           [~, ~, ext] = fileparts(reg_ver);
           if isempty(ext)
               reg_ver = [reg_ver, '.mat'];
           end
           file_name = sprintf('%s_%s_%s', dataset_name, stack, reg_ver);
           output = fullfile(obj.fp_metadata_folder(dataset_name, stack), file_name);            
        end
        
        function output = load_registration_data(obj, dataset_name, stack, reg_ver)
            output = obj.load_data(obj.fp_registration_data(dataset_name, stack, reg_ver));
        end


%% Graph
        function output = fp_graph_folder(obj, dataset_name, stack, create_dir)
            if nargin < 4
                create_dir = true;
            end
            output = fullfile(obj.fp_processed_data(dataset_name, stack), 'graph');
            if ~isfolder(output) && create_dir
                warning('Target folder does not exist. Create one. ');
                mkdir(output);
            end
        end
        
        function output = fp_graph_in_block_folder(obj, dataset_name, stack, grid_name, create_dir)
            if nargin < 5
                create_dir = true;
            end
            
            output = fullfile(obj.fp_graph_folder(dataset_name, stack), 'graph_in_block', grid_name);
            if ~isfolder(output) && create_dir
                warning('Target folder does not exist. Create one. ');
                mkdir(output);
            end
        end 
        
        function output = fp_graph_in_block_file(obj, dataset_name, stack, grid_name, idx_1, idx_2, layer, create_dir)
            if nargin < 8
                create_dir = true;
            end
            filename = sprintf('%s_%s_%s_graph_in_block_%d_%d_%d.mat', dataset_name, stack, grid_name, idx_1, idx_2, layer);
            output = fullfile(obj.fp_graph_in_block_folder(dataset_name, stack, grid_name, create_dir), filename);
        end
        
        function write_graph_in_block_v1(obj, block_segmentation)
            % ... some confusing notation...
            if isfield(block_segmentation, 'parameters')
                info_str = block_segmentation.parameters;
                layer = info_str.grid_layer;
                idx_1 = info_str.grid_idx_1;
                idx_2 = info_str.grid_idx_2;
            elseif isfield(block_segmentation, 'info')
                info_str = block_segmentation.info;
                if isfield(info_str, 'layer')
                    layer = info_str.layer;
                elseif isfield(info_str, 'grid_layer')
                    layer = info_str.grid_layer;
                end
                if isfield(info_str, 'idx_1')
                    idx_1 = info_str.idx_1;
                elseif isfield(info_str, 'grid_idx_1')
                    idx_1 = info_str.grid_idx_1;
                end
                if isfield(info_str, 'idx_2')
                    idx_2 = info_str.idx_2;
                elseif isfield(info_str, 'grid_idx_2')
                    idx_2 = info_str.grid_idx_2;
                end
            else
                error('Please specified block information in the input structure');
            end
            dataset_name = info_str.dataset_name;
            stack = info_str.stack;
            grid_name = info_str.grid_name;
            fp = obj.fp_graph_in_block_file(dataset_name, stack, grid_name, idx_1, idx_2, layer, true);
            save(fp, '-struct', 'block_segmentation', '-V7.3');
        end
        
        function write_graph_in_block(obj, block_segmentation, dataset_name, stack, graph_version, idx_1, ...
                idx_2, layer)
            fp = obj.fp_graph_in_block_file(dataset_name, stack, graph_version, idx_1, idx_2, layer, true);
            save(fp, '-struct', 'block_segmentation', '-V7.3');
        end
        
        
        function output = load_graph_in_block(obj, dataset_name, stack, grid_version, idx_1, idx_2, layer)
            fp = obj.fp_graph_in_block_file(dataset_name, stack, grid_version, idx_1, idx_2, layer, false);
            output = obj.load_data(fp);
        end
%%     Block data
        function output = fp_block_data_folder(obj, dataset_name, stack, version)
            output = fullfile(obj.fp_processed_data(dataset_name, stack),'block_data', version);            
        end
        
        function output = fp_block_data_file(obj, dataset_name, stack, version,idx1, idx2, layer)
            filename = sprintf('%s_%s_%s_block_data_%d_%d_%d.tiff', ...
                dataset_name, stack, version, idx1, idx2, layer);
            output = fullfile(obj.fp_block_data_folder(dataset_name, stack, version), filename);
        end
        
        function write_block_data_file(obj, data, dataset_name, stack, version, idx1, idx2, layer)
            filepath = obj.fp_block_data_file(dataset_name, stack, version, idx1, idx2, layer);
            [folder, ~, ~] = fileparts(filepath);
            if ~isfolder(folder)
                mkdir(folder);
            end
            obj.write_tiff_stack(data, filepath);
        end
        function output = load_block_data(obj, dataset_name, stack, version, idx1, idx2, layer)
            fp = obj.fp_block_data_file(dataset_name, stack, version, idx1, idx2, layer);
            output = obj.load_data(fp);   
        end
        
        function output = fp_block_data_DICOM_folder(obj, dataset_name, stack, version)
            output = fullfile(obj.fp_processed_data(dataset_name, stack), 'block_data_DICOM',version);
        end
        
        function output = fp_block_data_DICOM_file(obj, dataset_name, stack, version,idx1, idx2, layer)
            filename = sprintf('%s_%s_%s_block_data_%d_%d_%d.dcm', ...
                dataset_name, stack, version, idx1, idx2, layer);
            output = fullfile(obj.fp_block_data_DICOM_folder(dataset_name, stack, version), filename);
        end
        
        function write_block_data_DICOM_file(obj, data, dataset_name, stack, version, idx1, idx2, layer)
            % write_block_data_DICOM_file save 3d array DATA to designated
            % folder. The data is saved as [idx1, idx2, color, idx_frame]
            % in DICOM file. 
            filepath = obj.fp_block_data_DICOM_file(dataset_name, stack, version, idx1, idx2, layer);
            [folder, ~, ~] = fileparts(filepath);
            if ~isfolder(folder)
                warning('Target folder made.');
                mkdir(folder);
            end
            [lx, ly, lz] = size(data);
            dicomwrite(reshape(data, [lx, ly, 1, lz]), filepath);
        end
        
        function output = load_block_data_DICOM(obj, dataset_name, stack, version, idx1, idx2, layer)
            fp = obj.fp_block_data_DICOM_file(dataset_name, stack, version, idx1, idx2, layer);
            output = obj.load_data(fp);
        end
        
        function combined_block = load_blocks_data(obj, dataset_name, stack, grid_version, idx_X_list, idx_Y_list, layer_list, dtype)
            % One potential problem for this function is that the size of
            % the empty block might not be grid_info.block_size, if these
            % empty blocks are on the boundary of the whole data set. TO be
            % fixed later. 
            if nargin < 8
                dtype = 'uint16';
            end
            
            grid_info = obj.load_grid(dataset_name, stack, grid_version);
            num_layer = length(layer_list);
            num_X = length(idx_X_list);
            num_Y = length(idx_Y_list);

            block_cells = cell([num_X, num_Y, num_layer]);
            for layer_idx = 1 : num_layer
                layer = layer_list(layer_idx);
                for x_count = 1 : num_X
                    idx1 = idx_X_list(x_count);
                    for y_count = 1 : num_Y
                        idx2 = idx_Y_list(y_count);
                        if grid_info.bbox_xy_linear_idx_mat{layer}(idx1, idx2)
%                             fprintf('Loading data from %d %d %d\n', idx1, idx2, layer);
                            block_data = obj.load_block_data(dataset_name, stack, grid_version, idx1, idx2, layer);
                        else
                            warning('Block (%d, %d, %d) does not exist. Use zeros array.', idx1, idx2, layer);
                            if ~isfield(grid_info, 'grid2D')
                                block_data = zeros(grid_info.block_size, dtype);
                            else
                                block_data = zeros([grid_info.grid2D.ll(sub2ind(grid_info.grid2D.size, idx1,idx2),:), grid_info.gridZ.ll(layer)], dtype);
                            end
                        end
                        
                        if layer_idx ~= num_layer
                            block_data = block_data(:,:, 1:end-grid_info.block_overlap);
                        end
                        if x_count ~= num_X
                            block_data = block_data(1:end-grid_info.block_overlap, :, :);
                        end
                        if y_count ~= num_Y
                            block_data = block_data(:,1:end-grid_info.block_overlap,:);
                        end
                        block_cells{x_count, y_count, layer_idx} = block_data;
                    end
                end
            end
            combined_block = cell2mat(block_cells);
        end
        
        function output = load_block_data_by_bbox(obj, grid_info, bbox)
            % Input:
            %   grid_info: with field 'dataset_name', 'stack',
            %   'grid_version'(which specify the type of blocks to load),
            %   and grid size information. 
            %   bbox: [min_row, min_col, min_sec, row_lenght, col_length, sec_layer]
%             bbox = [3361, 1569, 1961, 232, 232, 64];
            min_row = bbox(1);
            min_col = bbox(2);
            min_sec = bow(3);
            row_length = bbox(4);
            col_length = bbox(5);
            sec_length = bbox(6);
            
            bbox_layer_idx0 = ceil(min_sec/(grid_info.block_z_size - grid_info.block_overlap));
            bbox_layer_idx1 = ceil( (min_sec + sec_length - 1 - grid_info.block_overlap)/(grid_info.block_z_size - grid_info.block_overlap));
            
            bbox_row_idx0 = ceil(min_row / (grid_info.block_xy_size - grid_info.block_overlap));
            bbox_row_idx1 = ceil( (min_row + row_length - 1 - grid_info.block_overlap) / (grid_info.block_xy_size - grid_info.block_overlap));
            
            bbox_col_idx0 = ceil(min_col / (grid_info.block_xy_size - grid_info.block_overlap));
            bbox_col_idx1 = ceil( (min_col + col_length - 1 - grid_info.block_overlap) / (grid_info.block_xy_size - grid_info.block_overlap));
            
            stitched_block_data = obj.load_blocks_data(grid_info.dataset_name, grid_info.stack, grid_info.grid_name, ...
                bbox_layer_idx0:bbox_layer_idx1, bbox_row_idx0:bbox_row_idx1, bbox_col_idx0:bbox_col_idx1);
            % Compute the local coordinate
            block_at_origin_idx_in_layer = grid_info.bbox_xy_linear_idx_mat{bbox_layer_idx0}(bbox_row_idx0, bbox_col_idx0);
            stitched_block_coordinate_origin = grid_info.bbox_xyz_mmll{bbox_layer_idx0}(block_at_origin_idx_in_layer,1:3) - 1;
            
            target_bbox_pos = bbox(1:3) - stitched_block_coordinate_origin;
            
            output = stitched_block_data(target_bbox_pos(1): target_bbox_pos(1) + row_length - 1, ...
                target_bbox_pos(2): target_bbox_pos(2) + col_length - 1, ...
                target_bbox_pos(3): target_bbox_pos(3) + sec_length - 1);
        
        end
        
        function output = fp_block_dt_folder(obj, dataset_name, stack, version)
            output = fullfile(obj.fp_processed_data(dataset_name, stack), 'block_dt',version);
        end
        
        function output = fp_block_dt_file(obj, dataset_name, stack, version, idx1, idx2, layer)
            filename = sprintf('%s_%s_%s_block_dt_%d_%d_%d.mat', ...
                dataset_name, stack, version, idx1, idx2, layer);
            output = fullfile(obj.fp_block_dt_folder(dataset_name, stack, version), filename);
        end
        
        function write_block_dt_file(obj, data, dataset_name, stack, version, idx1, idx2, layer)
            % write_block_data_DICOM_file save 3d array DATA to designated
            % folder. The data is saved as [idx1, idx2, color, idx_frame]
            % in DICOM file. 
            filepath = obj.fp_block_dt_file(dataset_name, stack, version, idx1, idx2, layer);
            [folder, ~, ~] = fileparts(filepath);
            if ~isfolder(folder)
                warning('Target folder made.');
                mkdir(folder);
            end
            save(filepath, '-struct', 'data');
        end
        
        function output = load_block_dt(obj, dataset_name, stack, version, idx1, idx2, layer)
            fp = obj.fp_block_dt_file(dataset_name, stack, version, idx1, idx2, layer);
            output = obj.load_data(fp);
        end
        
        % Skeleton 
        function output = fp_block_skl_folder(obj, dataset_name, stack, version)
            output = fullfile(obj.fp_processed_data(dataset_name, stack), 'block_skl',version);
        end
        
        function output = fp_block_skl_file(obj, dataset_name, stack, version, idx1, idx2, layer)
            filename = sprintf('%s_%s_%s_block_skl_%d_%d_%d.mat', ...
                dataset_name, stack, version, idx1, idx2, layer);
            output = fullfile(obj.fp_block_skl_folder(dataset_name, stack, version), filename);
        end
        
        function write_block_skl_file(obj, data, dataset_name, stack, version, idx1, idx2, layer)
            filepath = obj.fp_block_skl_file(dataset_name, stack, version, idx1, idx2, layer);
            [folder, ~, ~] = fileparts(filepath);
            if ~isfolder(folder)
                warning('Target folder made.');
                mkdir(folder);
            end
            save(filepath, '-struct', 'data');
        end
        
        function output = load_block_skl(obj, dataset_name, stack, version, idx1, idx2, layer)
            fp = obj.fp_block_skl_file(dataset_name, stack, version, idx1, idx2, layer);
            output = obj.load_data(fp);
        end        
        
        function output_str = load_blocks_skl_str(obj, dataset_name, stack, primary_version, idx_1_list, idx_2_list, layer_list, second_version)
            % One potential problem for this function is that the size of
            % the empty block might not be grid_info.block_size, if these
            % empty blocks are on the boundary of the whole data set. TO be
            % fixed later. 
            if nargin < 8
                second_version = [];
            end
            num_layer = length(layer_list);
            num_X = length(idx_1_list);
            num_Y = length(idx_2_list);
            block_cells = deal(cell([num_X, num_Y, num_layer]));
            cell_array_size = size(block_cells);
            num_blocks = num_layer * num_X * num_Y;
            block_global_bbox_mmxx_list = nan(num_blocks, 6);
            is_from_primary_Q = false(num_blocks, 1);
            for layer_idx = 1 : num_layer
                layer = layer_list(layer_idx);
                for x_count = 1 : num_X
                    idx1 = idx_1_list(x_count);
                    for y_count = 1 : num_Y
                        idx2 = idx_2_list(y_count);
                        
                        block_ind = sub2ind(cell_array_size, x_count, y_count, layer_idx);
                        file_fp = obj.fp_block_skl_file(dataset_name, stack, primary_version, idx1, idx2, layer);
                        if ~isfile(file_fp) && ~isempty(second_version)
                            fprintf('The primary version of the skeleton does not exist. Look for the secondary version.\n');
                            file_fp = obj.fp_block_skl_file(dataset_name, stack, second_version, idx1, idx2, layer);
                            if ~isfile(file_fp)
                                fprintf('The second version of the skeleton does not exist. Skip this block\n');
                            end
                        else
                            is_from_primary_Q(block_ind) = true;
                        end
                        tmp_str = obj.load_data(file_fp);
                        
                        block_cells{block_ind} = tmp_str;
                        block_global_bbox_mmxx_list(block_ind, :) = tmp_str.global_bbox_mmxx;
                    end
                end
            end
            array_gbbox_mmxx = nan(1, 6);
            array_gbbox_mmxx(1:3) = min(block_global_bbox_mmxx_list(:, 1:3), [], 1, 'omitnan');
            array_gbbox_mmxx(4:6) = max(block_global_bbox_mmxx_list(:, 4:6), [], 1, 'omitnan');
            array_bbox_ll = array_gbbox_mmxx(4:6) - array_gbbox_mmxx(1:3) + 1;
            block_bbox_mmxx_list = bsxfun(@minus, block_global_bbox_mmxx_list, ...
                [array_gbbox_mmxx(1:3), array_gbbox_mmxx(1:3)] - 1);            
            block_array_bbox_ll = block_global_bbox_mmxx_list(:, 4:6) - ...
                block_global_bbox_mmxx_list(:, 1:3) + 1;
            is_valid_bbox_Q = all(~isnan(block_global_bbox_mmxx_list), 2);
            [sub_cell, ind_cell, radius_cell, label_cell] = deal(cell(num_blocks, 1));            
            
            for iter_bbox = 1 : num_blocks
               % Convert indices in local 240 cubes to indices in the combined array
               if is_valid_bbox_Q(iter_bbox)
                   tmp_ind_240 = block_cells{iter_bbox}.ind;
                   tmp_sub_240 = fun_ind2sub(block_array_bbox_ll(iter_bbox, :), tmp_ind_240);
                   tmp_sub = bsxfun(@plus, tmp_sub_240, block_bbox_mmxx_list(iter_bbox, 1:3) - 1);
                   tmp_ind = sub2ind(array_bbox_ll, tmp_sub(:, 1), tmp_sub(:, 2), tmp_sub(:, 3));
                   ind_cell{iter_bbox} = tmp_ind;
                   sub_cell{iter_bbox} = tmp_sub;
                   radius_cell{iter_bbox} = block_cells{iter_bbox}.r;
                   if isfield(block_cells{iter_bbox}, 'label')
                       % The primary version may have the filed 'label',
                       % but not the secondary version. 
                       label_cell{iter_bbox} = block_cells{iter_bbox}.label;
                   end
               end
            end
            % Output structure
            output_str = struct;
            output_str.dataset_name = dataset_name;
            output_str.stack = stack;
            output_str.dataset_size = tmp_str.dataset_size;
            output_str.global_bbox_mmxx = array_gbbox_mmxx;
            output_str.global_bbox_mmll = [array_gbbox_mmxx(1:3), array_bbox_ll];
            output_str.block_size = array_bbox_ll;
            output_str.idx_1 = [];
            output_str.idx_2 = [];
            output_str.layer = [];
            output_str.info.primary_version = primary_version;
            output_str.info.secondary_version = second_version;
            % Combine primary first, then use the secondary 
            primary_ind = cat(1, ind_cell{is_from_primary_Q});
            primary_radius = cat(1, radius_cell{is_from_primary_Q});
            primary_label = cat(1, label_cell{is_from_primary_Q});
            % Remove the duplicated voxels in the overlapping region of the
            % bounding boxes
            [primary_ind, unique_idx, ~] = unique(primary_ind, 'stable');
            primary_radius = primary_radius(unique_idx);
            if ~isempty(primary_label)
                primary_label = primary_label(unique_idx);
            end
            
            is_from_secondary_Q = ~is_from_primary_Q & is_valid_bbox_Q;
            if any(is_from_secondary_Q)
                primary_bbox_mmxx = nan(1, 6);
                primary_bbox_mmxx(1:3) = min(block_bbox_mmxx_list(is_from_primary_Q, 1:3), [], 1);
                primary_bbox_mmxx(4:6) = max(block_bbox_mmxx_list(is_from_primary_Q, 4:6), [], 1);
                assert(all(~isnan(primary_bbox_mmxx)));
                
                secondary_sub = cat(1, sub_cell{is_from_secondary_Q});
                secondary_ind = cat(1, ind_cell{is_from_secondary_Q});
                secondary_radius = cat(1, radius_cell{is_from_secondary_Q});
                secondary_label = cat(1, label_cell{is_from_secondary_Q});
                
                [secondary_ind, unique_idx, ~] = unique(secondary_ind, 'stable');
                secondary_radius = secondary_radius(unique_idx);
                secondary_sub = secondary_sub(unique_idx, :);
%                 Only add the voxels that are not in the bounding box of
%                 the primary version
                is_not_in_primary_bbox_Q = ~fun_voxel_sub_in_bbox_mmxx_Q(secondary_sub, primary_bbox_mmxx);
                secondary_ind = secondary_ind(is_not_in_primary_bbox_Q);
                secondary_radius = secondary_radius(is_not_in_primary_bbox_Q);
                if ~isempty(secondary_label)
                    secondary_label = secondary_label(unique_idx);
                    secondary_label = secondary_label(is_not_in_primary_bbox_Q);
                elseif ~isempty(primary_label)
                    fprintf('The voxel label in the secondary version of skeleton is empty. Use 0 to pad the label vector\n');
                    secondary_label = zeros(size(secondary_ind), 'like', primary_label);
                else
                    secondary_label = [];
                end
                primary_ind = cat(1, primary_ind, secondary_ind);
                primary_radius = cat(1, primary_radius, secondary_radius);
                primary_label = cat(1, primary_label, secondary_label);
                % Extra information
                output_str.info.primary_bbox_mmxx = primary_bbox_mmxx;
                output_str.info.is_in_primary_bbox_Q = ~is_not_in_primary_bbox_Q;
            end
            output_str.ind = primary_ind;
            output_str.r = primary_radius;
            output_str.label = primary_label;
        end        
%%     ITK SNAP data
        function output = fp_itksnap_data_folder(obj, dataset_name, stack, what)
            output = fullfile(obj.fp_processed_data(dataset_name, stack),'itk', what);            
        end
        
        function output = fp_itksnap_file_string(obj, dataset_name, stack, version, idx1, idx2, layer)
            filename = sprintf('itk_%s_%s_%s_%d_%d_%d', ...
                dataset_name, stack, version, idx1, idx2, layer);
            output = fullfile(obj.fp_itksnap_data_folder(dataset_name, stack, version), filename);
        end
        
        
        function output = fp_itksnap_image_file(obj, dataset_name, stack, version, idx1, idx2, layer, compressedQ)
            if nargin < 8
                compressedQ = true;
            end
            if compressedQ
                output = sprintf('%s_image.nii.gz', obj.fp_itksnap_file_string(dataset_name, stack, version, idx1, idx2, layer));
            else
                output = sprintf('%s_image.nii', obj.fp_itksnap_file_string(dataset_name, stack, version, idx1, idx2, layer));
            end
        end
        
        function output = fp_itksnap_mask_file(obj, dataset_name, stack, version, idx1, idx2, layer, compressedQ)
            if nargin < 8
                compressedQ = true;
            end
            if compressedQ
                output = sprintf('%s_mask.nii.gz', obj.fp_itksnap_file_string(dataset_name, stack, version, idx1, idx2, layer));
            else
                output = sprintf('%s_mask.nii', obj.fp_itksnap_file_string(dataset_name, stack, version, idx1, idx2, layer)); 
            end
        end
        
        function output = load_itksnap_mask(obj, dataset_name, stack, version, idx1, idx2, layer)
            fp = obj.fp_itksnap_mask_file(dataset_name, stack, version, idx1, idx2, layer);
            output = obj.load_data(fp);
        end
        
        function output = load_itksnap_image(obj, dataset_name, stack, version, idx1, idx2, layer)
            fp = obj.fp_itksnap_image_file(dataset_name, stack, version, idx1, idx2, layer);
            output = obj.load_data(fp);
        end
        
        
        function write_itksnap_mask(obj, data, dataset_name, stack, version, idx1, idx2, layer, add_info, compressedQ)
            itk_fp_mask = obj.fp_itksnap_file_string(dataset_name, stack, version, idx1, idx2, layer);
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
                save_infoQ = false;
            end
            
            if ~isfile(annotation_fn)
                annotation_fn = fullfile(obj.fp_itksnap_data_folder(dataset_name, stack, 'annotated'), ...
                    annotation_fn);
                if ~isfile(annotation_fn)
                    error('%s does not exist...\n', annotation_fn);
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
            for layer_idx = 1 : grid_info.num_grid_layer
                for bbox_idx = 1 : length(grid_info.bbox_xyz_mmll{layer_idx})
                    if ~isfield(grid_info, 'downsample_rate')
                        cropped_cube = crop_bbox3(bbox_image, max(1,floor(grid_info.bbox_xyz_mmll{layer_idx}(bbox_idx,:)/downsample_stack_info.downsample_ratio)),'default');
                    else
                        cropped_cube = crop_bbox3(bbox_image, max(1,floor(grid_info.bbox_xyz_mmll{layer_idx}(bbox_idx,:) * grid_info.downsample_rate /downsample_stack_info.downsample_ratio)),'default');
                    end
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
                    save_fn = fullfile(fp, sprintf('%s_block_info_for_%s.mat', fn, grid_info.grid_name));
                    blocks_info.filepath = save_fn;
                    fprintf('Block information saved to %s\n', save_fn);
                    save(save_fn, '-struct', 'blocks_info');
                end
            else
                disp('No blocks found, please check the input information...');
            end            
        end
        
        function blocks_info = load_annotation_blocks_info(obj, dataset_name, stack, fp)
            if ~isfile(fp)
                fp = fullfile(obj.fp_itksnap_data_folder(dataset_name, stack, 'annotated'), fp);
                if ~isfile(fp)
                    error('%s does not exist...\n', fp);
                end                
            end
            blocks_info = obj.load_data(fp);
            if isfield(blocks_info, 'blocks_info')
                blocks_info = blocks_info.blocks_info;
            end
        end
        
%%     Synchronize between server and the local workstation
     function save_to_server_and_delete_in_root(obj, fp,  no_wait)
         if nargin < 3
             no_wait = true;
         end
            
         [folder_path, ~, ~] = fileparts(fp);
         target_file_folder = strrep(folder_path, obj.ROOT_PATH, obj.SERVER_ROOT_PATH);
         if ~isfolder(target_file_folder)
             warning('Target folder does not exist. Create one');
             mkdir(target_file_folder)
         end
         if no_wait 
             cmd_str = sprintf('rsync -a --remove-source-files %s %s&', fp, target_file_folder);
         else
             cmd_str = sprintf('rsync -a --remove-source-files %s %s', fp, target_file_folder);
         end
         system(cmd_str);
     end
     function save_to_server_and_delete_in_scratch(obj, fp,  no_wait)
         if nargin < 3
             no_wait = true;
         end
            
         [folder_path, ~, ~] = fileparts(fp);
         target_file_folder = strrep(folder_path, obj.SCRATCH_ROOT_PATH, obj.SERVER_ROOT_PATH);
         if ~isfolder(target_file_folder)
             warning('Target folder does not exist. Create one');
             mkdir(target_file_folder)
         end
         if no_wait 
             cmd_str = sprintf('rsync -a --remove-source-files %s %s&', fp, target_file_folder);
         else
             cmd_str = sprintf('rsync -a --remove-source-files %s %s', fp, target_file_folder);
         end
         system(cmd_str);
     end

     function save_to_server_from_root(obj, fp, no_wait)
         if nargin < 3
             no_wait = true;
         end

         [folder_path, ~, ~] = fileparts(fp);
         target_file_folder = strrep(folder_path, obj.ROOT_PATH, obj.SERVER_ROOT_PATH);
         if ~isfolder(target_file_folder)
             warning('Target folder does not exist. Create one');
             mkdir(target_file_folder)
         end
         if no_wait 
             cmd_str = sprintf('rsync -a %s %s&', fp, target_file_folder);
         else
             cmd_str = sprintf('rsync -a %s %s', fp, target_file_folder);
         end
         system(cmd_str);
     end
     
     function save_to_server_from_stratch(obj, fp, no_wait)
         if nargin < 3
             no_wait = true;
         end

         [folder_path, ~, ~] = fileparts(fp);
         target_file_folder = strrep(folder_path, obj.SCRATCH_ROOT_PATH, obj.SERVER_ROOT_PATH);
         if ~isfolder(target_file_folder)
             warning('Target folder does not exist. Create one');
             mkdir(target_file_folder)
         end
         if no_wait 
             cmd_str = sprintf('rsync -a %s %s&', fp, target_file_folder);
         else
             cmd_str = sprintf('rsync -a %s %s', fp, target_file_folder);
         end
         system(cmd_str);
     end
    end
%% Methods for getting metadata
    methods(Static)

    end
    
    methods
        function output = n_raw_data_image_stack_size(obj, dataset_name, stack)
            output = [obj.n_raw_data_image_size(dataset_name, stack), obj.n_raw_data_image(dataset_name, stack) ];
        end
        
        function output = list_raw_data_image_section(obj, dataset_name, stack)
            output = 0 : obj.n_raw_data_image(dataset_name, stack)-1;
        end

        function output = invalid_section(obj, dataset_name, stack, section)
            output = ismember(section, obj.invalid_section_list(dataset_name, stack));
        end
        
        function output = valid_section_list(obj, dataset_name, stack)
            section_list = 1 : obj.n_raw_data_image(dataset_name, stack);
            output = setdiff(section_list, obj.invalid_section_list(dataset_name, stack));
        end
        
        function output = fp_dataset_info(obj, dataset_name, stack)
           output = fullfile(obj.fp_processed_data(dataset_name, stack), 'data_info.mat');
        end
        
        
        function output = load_dataset_info(obj, dataset_name, stack)
            fp = obj.fp_dataset_info(dataset_name, stack);
            output = obj.load_data(fp);
        end
        
        function write_dataset_info(obj, data_info)
            dataset_name = data_info.dataset_name;
            stack = data_info.stack;
            save(obj.fp_dataset_info(obj, dataset_name, stack), '-struct', 'data_info');
        end
        
    end   
%% Methods for loading files
    methods
        function output = load_data(obj, fp)
            if obj.opt_load_from_scratch_first_Q
               fp_scratch = strrep(fp, obj.ROOT_PATH, obj.SCRATCH_ROOT_PATH);
               if isfile(fp_scratch)
                   fp = fp_scratch;
               else
                   fprintf('Target file does not exist in the scratch folder\n');
               end
            end            
            fp = obj.fp_file_exist_locally(fp);
            [~, ~, ext] = fileparts(fp);
            switch ext
                case '.mat'
                    output = load(fp);
                    field_name = fieldnames(output);
                    if numel(field_name) == 1
                        output = getfield(output, field_name{1}); %#ok<GFLD>
                    end
                    
                case {'.tiff', '.tif'}
                    output = obj.load_single_tiff(fp);
                case {'.nii', '.gz'}
                    if strcmp(ext, '.gz')
                        assert(endsWith(fp, '.nii.gz'));
                    end
                    output = niftiread(fp);
                case {'.xml'}
                    output = xml_read(fp);
                case {'.nrrd'}
                    output = nrrdread(fp);
                case {'dcm'}
                    output = dicomread(fp);
                otherwise
                    error('Cannot recoginize file format %s', ext);
            end
        end
        
        function [target_folder, new_file_path_cell_array] = copy_to_scratch(obj, file_path_cell_array)
            target_folder = fullfile(obj.SCRATCH_ROOT_PATH, tempname);
            new_file_path_cell_array = cell(numel(file_path_cell_array),1);
            mkdir(target_folder);
            for fp_idx = 1 : numel(file_path_cell_array)
                [~, tmpFn, tmpExt] = fileparts(file_path_cell_array{fp_idx});
                new_file_path_cell_array{fp_idx} = fullfile(target_folder, [tmpFn, tmpExt]);
                copyfile(file_path_cell_array{fp_idx}, target_folder);
            end
        end
    end
%% Method for calling external library
%%      ITKSNAP
    methods
        function output = fp_itksnap_exe(obj)
            switch obj.computer_name
                case 'MACHINE 1'
                    output = fullfile(obj.EXTERNAL_LIB_ROOT_PATH, 'itksnap-3.8.0','bin','itksnap');
                case 'MACHINE 2'
                    output = '"C:\Program Files\ITK-SNAP 3.6\bin\ITK-SNAP.exe"';
            end 
            
        end
        

        function itk_fp_mask = visualize_itksnap(obj, data, mask, itk_file_string, saveImageQ)
            % Input:
            %   data: image stack, numerical array
            %   mask: logical array. If it is logical, convert it into
            %   integer, since niftiwrite doesn't support output binary array
            %   fp: output filepath, extension doesn't matter. If not
            %   given, saved to 'tmp' folder. 
            %   saveImageQ:
            % Output:
            %   Call ITKSNAP
            %   Return filepath for the mask 
            if nargin < 2
                itk_path = obj.fp_itksnap_exe();
                str_cmd = sprintf('%s &', itk_path);
                system(str_cmd);
                return;
            end
                
            if nargin < 4
                file_name = datestr(now,'mm_dd_yy_hh_MM_SS');
%                 folder_name_data = '/scratch';
%                 folder_name_mask = '/scratch';
                folder_name_data = obj.fp_temporary_folder();
                folder_name_mask = obj.fp_temporary_folder();
            elseif nargin <5
                [folder_name_mask, file_name, ~] = fileparts(itk_file_string);
                folder_name_data = obj.fp_temporary_folder();
            elseif saveImageQ
                [folder_name_mask, file_name, ~] = fileparts(itk_file_string);
                folder_name_data = folder_name_mask;
            else
                [folder_name_mask, file_name, ~] = fileparts(itk_file_string);
                folder_name_data = obj.fp_temporary_folder();
            end
            
            if ~isfolder(folder_name_mask)
                warning('Target folder made');
                mkdir(folder_name_mask);
            end
            
            if isnumeric(data)
                itk_fp_data = fullfile(folder_name_data, sprintf('%s_image.nii', file_name));
                if ~isfile(itk_fp_data)
                    niftiwrite(data, itk_fp_data);
                else
                fprintf('Use the existing image data to save time...\n');
                end
            elseif ischar(data)
                itk_fp_data = data;
                if ~isfile(itk_fp_data)
                    error('Image %s does not exist\n', itk_fp_data);
                end
            else
                error('data should be either numerical array or filepath');
            end
            
            if ~ischar(mask)
                itk_fp_mask = fullfile(folder_name_mask, sprintf('%s_mask', file_name));
                 if ~islogical(mask)
                    mask_max = max(mask(:));
                    mask_min = min(mask(:));
                    if (0 <= mask_min) && (mask_min <= mask_max) && (mask_max <= intmax('uint8'))
                        mask = uint8(mask);
                    elseif (0 <= mask_min) && (mask_min <= mask_max) && (mask_max <= intmax('uint16'))
                        mask = uint16(mask);
                    elseif (0 <= mask_min) && (mask_min <= mask_max) && (mask_max <= intmax('uint32'))
                        mask = uint32(mask);
                    else
                        mask = uint64(mask);
                    end
                else
                    mask = uint8(mask);
                 end
                 niftiwrite(mask, itk_fp_mask, 'Compressed', true);
            else
                itk_fp_mask = mask;
                if isfile(itk_fp_mask)
                    error('mask should be either numerical/logical array or filepath');
                end
            end
            itk_path = obj.fp_itksnap_exe();
            str_cmd = sprintf('%s -g %s -s %s.nii.gz &', itk_path, itk_fp_data, itk_fp_mask);
            system(str_cmd);
        end
        
        function visualize_block_location_in_itksnap(obj, dataset_name, stack, mask_version, grid_version, idx1, idx2, layer,downsample_ratio)
            if nargin < 9
                downsample_ratio = 8;
            end
            grid_info = obj.load_grid(dataset_name, stack, grid_version);
%             mask_info = obj.load_mask_info(dataset_name, stack, mask_version);
            downsampled_stack_info = obj.load_downsampled_stack_info(dataset_name, stack, mask_version, downsample_ratio);
            disp('Determine block location');
            tmp_match_1 = grid_info.bbox_xy_mm_id_pos{layer}(:,1) == idx1;
            tmp_match_2 = grid_info.bbox_xy_mm_id_pos{layer}(:,2) == idx2;
            block_linear_idx = find(and(tmp_match_1, tmp_match_2));
            if ~isfield(grid_info, 'downsample_rate')
                block_bbox_xy = round(grid_info.bbox_xy_mmxx{layer}(block_linear_idx,:)/downsample_ratio); 
                block_bbox_z = round(grid_info.bbox_z_mmxx{layer}/downsample_ratio);
            elseif isscalar(grid_info.downsample_rate)
                block_bbox_xy = round(grid_info.downsample_rate * grid_info.bbox_xy_mmxx{layer}(block_linear_idx,:)/downsample_ratio); 
                block_bbox_z = round(grid_info.downsample_rate * grid_info.bbox_z_mmxx{layer}/downsample_ratio);
            else
                error('Todo: anisotropic downsample rate');
            end
            
            disp('Construct block mask');
            vis_block_pos = zeros(downsampled_stack_info.stack_size, 'uint8');
            vis_block_pos(block_bbox_xy(1):block_bbox_xy(3), block_bbox_xy(2):block_bbox_xy(4), block_bbox_z(1):block_bbox_z(2)) = 1;
            itk_fp_mask = fullfile(obj.fp_temporary_folder(), sprintf('tmp_block_mask_%s', datestr(now,'mm_dd_yy_hh_MM_SS')));
            fprintf('Write itk mask to temporary folder: %s\n', itk_fp_mask);
            niftiwrite(vis_block_pos, itk_fp_mask, 'Compressed', true);
            if downsampled_stack_info.quick_access
                itk_fp_data = downsampled_stack_info.quick_access_nii_fp;
            else
                error('To be implemented');
                % TO DO
                % For the nii file not stored in the quick_access as nii
                % file
            end
            itk_path = obj.fp_itksnap_exe();
            str_cmd = sprintf('%s -g %s -s %s.nii.gz &', itk_path, itk_fp_data, itk_fp_mask);
            system(str_cmd); 
        end
        
    end
%% Method for output visualization result
    methods
        function fp = fp_analysis_image(obj, fn, subfolder_name)
            if nargin < 3
                subfolder_name = [];
            end
            
            if isempty(subfolder_name)
                fp = fullfile(obj.AnalysisResultRootPath, 'Figures', fn);
            else
                folder_path = fullfile(obj.AnalysisResultRootPath, 'Figures', subfolder_name);
                if ~isfolder(folder_path)
                    warning('Folder does not exist. Create folder');
                    mkdir(folder_path);
                end
                fp = fullfile(folder_path, fn);
            end               
               
        end
        
        function output = fp_visualization_folder(obj, dataset_name, stack)
            output = fullfile(obj.fp_Dataset(dataset_name), stack,'visualization');
        end
        
        function output = fp_visualization_data_file(obj, dataset_name, stack, what, ext)
            if nargin < 6
                ext = 'mat';
            end
            output = fullfile(obj.fp_visualization_folder(dataset_name, stack), ...
                sprintf('%s_%s_%s.%s', dataset_name, stack, what, ext));
        end
        
        function output = write_visualization_data(obj, data, dataset_name, stack, what, ext)
            if nargin < 6
                ext = 'mat';
            end
            fp = obj.fp_visualization_data_file(dataset_name, stack, what);
            [folder_path, ~, fp_ext]= fileparts(fp);
            if ~isempty(fp_ext)
                fp = sprintf('%s.%s', fp, ext);
            else
                ext = fp_ext;
            end
            if ~isfolder(folder_path)
               warning('Target folder does not exist. Create one');
               mkdir(folder_path);
            end
            
            switch ext
                case 'mat'
                    save(fp, '-struct', 'data', '-v7.3');
                otherwise 
                    error('Does not support. To be implemented');
            end
            output = 0;
        end
        
        function data = load_visualization_data(obj, dataset_name, stack, what, ext)
            if nargin < 5
                ext = 'mat';
            end
            fp = obj.fp_visualization_data_file(dataset_name, stack, what);
            [~, ~, fp_ext] = fileparts(fp);
            if isempty(fp_ext)
                fp = sprintf('%s.%s', fp, ext);
            else
                ext = fp_ext;
            end
            switch ext
                case 'mat'
                    data = obj.load_data(fp);
                otherwise
                    error('To be implemented');
            end
        end
        
        function fp = fp_visualization_image(obj, dataset_name, stack, fn, subfolder_name)
            if nargin < 5
                subfolder_name = [];
            end
            
            if isempty(subfolder_name)
                folder_path = obj.fp_visualization_folder(dataset_name, stack);
            else
                folder_path = fullfile(obj.fp_visualization_folder(dataset_name, stack), subfolder_name);
                if ~isfolder(folder_path)
                    warning('Folder does not exist. Create folder');
                    mkdir(folder_path);
                end
            end
            fp = fullfile(folder_path, fn);
        end        
    end
%% Read/write file 
    methods(Static)
        function clear_scratch_tmp_folder(target_folder, nowait)
            if nargin < 2
                nowait = true;
            end
            
            if contains(target_folder, 'scratch')
                if nowait
                   system(sprintf('rm -r %s&', target_folder)); 
                else
                   system(sprintf('rm -r %s', target_folder));
                end
            else
                error('Folder does not exist in the scratch folder');
            end
        end
        
        function output = load_single_tiff(data_fp, section_list, no_warningQ)
            warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
            warning('off', 'imageio:tiffmexutils:libtiffWarning');
            % load_single_tiff loads tiff (stack) at data_fp. If
            % section_list is not given, loads all section in the stack. 
            if nargin < 2
                section_list = nan;
                no_warningQ = false;
            elseif nargin < 3
                no_warningQ = false;
            end
            
            try
                n_specified_slice = length(section_list);
                image_info = imfinfo(data_fp);
                tifLink = Tiff(data_fp, 'r');
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
                    case {8, 24}
                        image_type = 'uint8';
                    case 1
                        image_type = 'logical'; % In matlab, logical is actually uint8
                        if ~ no_warningQ
                            warning('Read logical array through Tiff library, which inverts the value in image output by imwrite. Check if result is consistant with your expectation.');
                        end
                    otherwise
                        error('Unrecongnized image bit depth %d', image_bit_depth);
                end
                
                if isnan(section_list)
                    if image_bit_depth ~= 24
                        output = zeros(image_size_1, image_size_2, num_slice, image_type);
                        for iSection  = 1 : num_slice
                            tifLink.setDirectory(iSection);
                            output(:,:,iSection) = tifLink.read();
                        end
                    else
                        output = zeros(image_size_1, image_size_2, 3, num_slice, image_type);
                        for iSection = 1 : num_slice
                            tifLink.setDirectory(iSection);
                            output(:,:,:,iSection) = tifLink.read();
                        end
                    end

                elseif n_specified_slice == 1
                    output = tifLink.read();
                else
                    if image_bit_depth ~= 24
                        output = zeros(image_size_1, image_size_2, n_specified_slice, image_type);
                        for iSection = 1 : n_specified_slice
                            tmp_read_sec = section_list(iSection);
                            tifLink.setDirectory(tmp_read_sec);
                            output(:,:,iSection) = tifLink.read();                           
                        end
                    else
                        output = zeros(image_size_1, image_size_2, 3, n_specified_slice, image_type);
                        for iSection = 1 : n_specified_slice
                            tmp_read_sec = section_list(iSection);
                            tifLink.setDirectory(tmp_read_sec);
                            output(:,:,iSection) = tifLink.read();                           
                        end
                    end
                end
            catch ME
                error(ME.message)
            end 
            tifLink.close();
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
        
        function write_tiff_stack(inputArray, fp, image_type,output_mode)
            % Adapt from YoonOh Tak's scrip saveastiff
            % This function also support output single section tiff.
            [folder_path, ~, ~] = fileparts(fp);
            if ~isfolder(folder_path)
                warning('Folder does not exist. Create folder.');
                mkdir(folder_path)
            end
            if nargin < 3
                image_type = 'grayscale';
                output_mode = 'overwrite';
            elseif nargin < 4
                output_mode = 'overwrite';
            end
            
            switch image_type
                case 'grayscale'
                    tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                    [im_height, im_width, im_depth] = size(inputArray);
                    tagstruct.SamplesPerPixel = 1;
                case {'color', 'RGB'}
                    tagstruct.Photometric = Tiff.Photometric.RGB;
                    [im_height, im_width, im_cc, im_depth] = size(inputArray);
                    tagstruct.SamplesPerPixel = im_cc;
                    if im_cc == 4
                        tagstruct.ExtraSamples = Tiff.ExtraSamples.AssociatedAlpha;
                    end
            end
            tagstruct.ImageLength = im_height;
            tagstruct.ImageWidth = im_width;
            tagstruct.Compression = Tiff.Compression.None;
            % (RGB RGB,RGB RGB,RGB RGB),
            % http://www.awaresystems.be/imaging/tiff/tifftags/planarconfiguration.html
            tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
            switch class(inputArray)
                case {'uint8', 'uint16', 'uint32', 'logical'}
                    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
                case {'int8', 'int16', 'int32'}
                    tagstruct.SampleFormat = Tiff.SampleFormat.Int;
                case {'single', 'double', 'uint64', 'int64'}
                    tagstruct.SampleFormat = Tiff.SampleFormat.IEEEFP;
                otherwise
                    error('Unsupported format');
            end
            % Bits per sample
            switch class(inputArray)
                case {'logical'}
                    tagstruct.BitsPerSample = 1;
                case {'uint8', 'int8'}
                    tagstruct.BitsPerSample = 8;
                case {'uint16', 'int16'}
                    tagstruct.BitsPerSample = 16;
                case {'uint32', 'int32'}
                    tagstruct.BitsPerSample = 32;
                case {'single'}
                    tagstruct.BitsPerSample = 32;
                case {'double', 'uint64', 'int64'}
                    tagstruct.BitsPerSample = 64;
                otherwise
                    error('Unsupported format');
            end
            % Rows per strip
            maxstripsize = 8192;
            byte_depth = tagstruct.BitsPerSample/8;
            tagstruct.RowsPerStrip = ceil(maxstripsize/(im_width*(byte_depth)*tagstruct.SamplesPerPixel)); % http://www.awaresystems.be/imaging/tiff/tifftags/rowsperstrip.html
            inputArray = reshape(inputArray, im_height, im_width, tagstruct.SamplesPerPixel, im_depth);
            % Open file
            switch output_mode
                case {'overwrite', 'w'}
                    if (im_height * im_width * tagstruct.SamplesPerPixel * im_depth * byte_depth) > 2^32 - 1
                        tfile = Tiff(fp, 'w8');
                    else
                        tfile = Tiff(fp, 'w');
                    end
                case {'append', 'r+'}
                    tfile = Tiff(fp, 'r+');
                    while ~tfile.lastDirectory() % Append a new image to the last directory of an exiting file
                        tfile.nextDirectory();
                    end
                    tfile.writeDirectory();
            end
            for secID = 1 : im_depth
                tfile.setTag(tagstruct);
                tfile.write(inputArray(:,:,:,secID));
                if secID ~= im_depth
                    tfile.writeDirectory();
                end
            end
            tfile.close();
        end
        
        function write_image_stack_as_avi(image_stack, fp)
            [folder, ~, ~] = fileparts(fp);
            if ~isfolder(folder)
                warning('Target folder does not exist. Create %s', folder);
                mkdir(folder);
            end
            avi_str = VideoWriter(fp);
            avi_str.FrameRate = 10;
            avi_str.Quality = 75;
            open(avi_str)
            for frame_idx = 1 : size(image_stack, 3)
                writeVideo(avi_str, image_stack(:,:,frame_idx));
            end
        end    
        
        function write_data(output_fp, data)
           [tmp_f, tmp_n, tmp_ext] = fileparts(output_fp);
           
           if ~isfolder(tmp_f)
               mkdir(tmp_f);
           end
           switch tmp_ext
               case {'.mat'}
                   if isstruct(data)
                       save(output_fp, '-struct', 'data', '-v7.3');
                   else
                       save(output_fp, 'data', '-v7.3');
                   end
               case {'.tiff', '.tif'}
                   write_tiff_stack(data, output_fp);
               otherwise
                   error('Unrecognized output file format');
           end
           assert(isfile(output_fp), sprintf('Fail to write file %s for unknown reason', ...
               output_fp));
        end
    end
%% Program Metadata
    methods
        function output = fp_constructor_for_script(obj, subfolder_name)
            output = fullfile(obj.SCRIPT_PATH, subfolder_name);
        end
        
        % Filepath to the parameters XML file for tasks
        function output = fp_task_parameters(obj, varargin)
            output = fullfile(obj.SCRIPT_PATH, 'Task', 'task_parameters', varargin{:});
        end
        
        function output = load_task_parameters(obj, varargin)
            fp = obj.fp_task_parameters(varargin{:});
            output = obj.load_data(fp);
        end
        
        function exit_code = write_task_parameters(obj, data, varargin)
            assert(isstruct(data), 'The data to be saved should be a MATLAB structure');
            data_fp = obj.fp_task_parameters(varargin{:});
            [folder, ~, ext] = fileparts(data_fp);
            if ~isfolder(folder)
                mkdir(folder);
            end
            
            switch ext
                case '.xml'
                    xml_write(data_fp, data);
                case '.mat'
                    save(data_fp, '-struct', 'data');
                otherwise
                    error('Unrecognized file type');
            end
            exit_code = 0;
            fprintf('Write data to %s\n', data_fp);
        end
    end
%% Parallel computing metadata
    methods
        function fp = fp_task_root_folder(obj, dataset_name, stack)
            fp = fullfile(obj.fp_Dataset(dataset_name), stack, 'task');
        end
        
        function fp = fp_task_folder(obj, dataset_name, stack, task_name)
            fp = fullfile(obj.fp_task_root_folder(dataset_name, stack), task_name);
        end
        
        function fp = fp_task_file(obj, dataset_name, stack, task_folder, task_name)
            fp = fullfile(obj.fp_task_folder(dataset_name, stack, task_folder), task_name);
            folder_name = fileparts(fp);
            if ~isfolder(folder_name)
                mkdir(folder_name);
            end
        end    
                % Write mat file to Scratch folder and its subfolder
        function file_name = write_mat_file_to_scratch_folder(obj, file_name, data)
            % Check if the file name is an absolute path 
            if ~startsWith(file_name, filesep)
                file_name = fullfile(obj.SCRATCH_ROOT_PATH, 'tmp', file_name);
            end
            [target_folder, ~, target_ext] = fileparts(file_name);
            if ~isfolder(target_folder)
                mkdir(target_folder);
            end
            switch target_ext
                case {'.mat'}
                    if isstruct(data)
                        save(file_name, '-struct', 'data');
                    else
                        save(file_name, 'data');
                    end
                otherwise 
                    error('Unrecognized output data type');
            end
        end
    end
%% Dataset specific methods - MouseLight
    methods
        function output = mouselight_load_octree_data(obj, stack, octree_coordinate)
            raw_data_folder = obj.fp_raw_data_folder('WholeBrain', stack);
            if nargin < 3 || isempty(octree_coordinate)
                file_path = fullfile(raw_data_folder, 'default.0.tif');
            else
                file_path = fullfile(raw_data_folder, join(string(octree_coordinate), filesep), 'default.0.tif');
            end
            if ~isfile(file_path)
                file_path = obj.fp_file_exist_locally(file_path);
                if ~isfile(file_path)
                    error('File does not exist');
                end
            end
            warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
            output = obj.load_data(file_path);
            warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffWarning');
        end
        
    end
%% Allen Atlas
    methods
        function data_str = load_single_region_stat_data_by_atlas_id(obj, dataset_name, stack_list, region_atlas_id)
            % The 'Allen_atlas.mat' is in the Registration folder of the code.
            persistent registration_str;
            if isempty(registration_str)
                registration_str = obj.load_data('Allen_atlas.mat');
            end            
            region_list_ind = full(registration_str.id_2_ind(region_atlas_id));
            region_name = registration_str.structure_table.name(region_list_ind);
            region_name = region_name{1};
            region_name_abv = registration_str.structure_table.acronym(region_list_ind);
            region_name_abv = region_name_abv{1};
            num_stack = numel(stack_list);
            stack_data = cell(num_stack, 1);
            valid_data_Q = false(num_stack, 1);
            for iter_stack = 1 : num_stack
                tmp_stack = stack_list{iter_stack};
                region_data_folder = fullfile(obj.fp_metadata_folder(dataset_name, ...
                    tmp_stack), 'region_data');
                
                region_data_filepath = fullfile(region_data_folder, ...
                    sprintf('%s_%s_region_data_%s.mat', dataset_name, tmp_stack, ...
                    strrep(region_name, ' ', '_')));
                
                if isfile(region_data_filepath)
                    tmp_tic = tic;
                    stack_data{iter_stack} = obj.load_data(region_data_filepath);
                    fprintf('Finish loading %s. Elapsed time is %f seconds\n', ...
                        region_data_filepath, toc(tmp_tic));
                    if isfield(stack_data{iter_stack}, 'cc_1') && isfield(stack_data{iter_stack}, 'cc_2')
                        valid_data_Q(iter_stack) = true;
                    end
                else
                    fprintf('File %s does not exist\n', region_data_filepath);
                end
            end
            data_str = struct;
            data_str.dataset_name = dataset_name;
            data_str.stack_list = stack_list;
            data_str.region_atlas_id = region_atlas_id;
            data_str.region_atlas_ind = region_list_ind;
            data_str.region_name = region_name;
            data_str.region_name_abv = region_name_abv;
            data_str.num_stack = num_stack;
            data_str.is_valid_data_Q = valid_data_Q;
            data_str.stack_data = stack_data;
        end      
    end
%% Parallel computing control
    methods(Static)
        function run_command_on_machine(machine_name, cmd_str, no_wait_Q)
            if nargin < 3
                no_wait_Q = true;
            end
            % Use the following information to setup the machine first:
            % 1. https://www.ostechnix.com/how-to-create-ssh-alias-in-linux/
            % 2. https://askubuntu.com/questions/8653/how-to-keep-processes-running-after-ending-ssh-session
            ssh_cmd_str = strjoin({'"', cmd_str, '"'}, ' ');
            if no_wait_Q
                ssh_str = strjoin({'ssh', machine_name, ssh_cmd_str, '&'}, ' ');
            else
                ssh_str = strjoin({'ssh', machine_name, ssh_cmd_str}, ' ');
            end
            system(ssh_str);
        end
    end
    
    methods 
        function sync_script_to_server(obj)
            sync_source = obj.SCRIPT_PATH;
            sync_target = obj.SERVER_SCRIPT_PATH;
            % Do not sync the hidden files (e.g. under .git/)
            sync_str = sprintf('rsync -rav --exclude=".*" %s/ %s', ...
                sync_source, sync_target);
            system(sync_str);
            fprintf('Finish synchronizing local script to the data folder.\n');
        end
    end
end