%% Load the Allen Atlas
set_env;
DataManager = FileManager;
atlas_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_25.nrrd');
annotation_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd');
atlas_im = DataManager.load_data(atlas_fp);
atlas_annotation = DataManager.load_data(annotation_fp);
DataManager.visualize_itksnap(atlas_im, atlas_annotation);
allen_atlas = load('Allen_atlas.mat');
%% Load registered data
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124'};
num_stack = numel(stack_list);
registration_version = 'Allen_2017_25um_nonrigid.mat';
registration_cell = cell(num_stack, 1);
for iter_stack = 1 : num_stack
    registration_cell{iter_stack} = DataManager.load_registration_data(dataset_name, ...
        stack_list{iter_stack}, registration_version);
end
%%
analysis_region_id_list = {453 961 1080 477 549 [302 294] 4 771 354 528 343, 776}; % Remove hypothalamus, use hippocampus region
num_region = numel(analysis_region_id_list);
analysis_region_name = cell(num_region, 1);
for iter_region = 1 : num_region
    tmp_atlas_id = analysis_region_id_list{iter_region};
    if isscalar(tmp_atlas_id)
        analysis_region_name{iter_region} = allen_atlas.structure_table.name{full(allen_atlas.id_2_ind(tmp_atlas_id))};
    end
end
analysis_region_name{6} = 'Superior colliculus';
%%
tmp_region_idx = 1;
tmp_region_id = analysis_region_id_list{tmp_region_idx};
[subregion_cc, subregion_mask] = fun_registration_get_region_cc(...
    registration_cell{iter_stack}, tmp_region_id);
%% Allen atlas volume
% subregion_list_ind = 