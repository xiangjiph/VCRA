% Load the whole brain mask and slipt it into pieces 
wb_im = DataManager.load_single_tiff('/data/Vessel/WholeBrain/ML_2018_08_15/processed_data/whole_stack_d16x.tiff');
wb_im = permute(wb_im, [1, 3, 2]);
tmp_folder = fullfile('/data/Vessel/WholeBrain/ML_2018_08_15/processed_data', 'whole_stack_d16x_sagittal');
if ~isfolder(tmp_folder)
    mkdir(tmp_folder);
end
for iter_layer = 1 : size(wb_im, 3)
    fprintf('Writing tiff stack %d\n', iter_layer);
    DataManager.write_tiff_stack(wb_im(:, :, iter_layer), fullfile(tmp_folder, sprintf('whole_stack_d16x_sagittal_section_%d.tiff', iter_layer)));
end
