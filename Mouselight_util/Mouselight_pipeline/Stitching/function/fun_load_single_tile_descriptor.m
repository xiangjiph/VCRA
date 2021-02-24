function des = fun_load_single_tile_descriptor(descriptorfolder, relative_path, channel, ext)
% fun_load_single_tile_descriptor loads the descriptor of the tile. It is
% part of fun_xymatch_vessel_with_IO, which load the descriptor of the tile
% on the fly, instead of loading all the descriptor before calling the xy
% matching function. 
% 
% Note: 
% 1. Support for multiple channel descriptor to be implemented. 
% 2. Descriptor pre-selection can be added later. 
%
%
if nargin < 4
    ext = '.mat';
end
if iscell(channel)
    tmp_descriptor_name = sprintf('*.%s%s', channel{1}, ext);
elseif isnumeric(channel)
     tmp_descriptor_name = sprintf('*.%d%s', channel, ext);
elseif isstring(chennel)
    tmp_descriptor_name = sprintf('*.%s%s', channel, ext);
end
des_center_fp_str = dir(fullfile(descriptorfolder, relative_path, tmp_descriptor_name));
if isempty(des_center_fp_str.name)
    fprintf('Missing descriptor for tile %s\n', relative_path);
    des = [];
else
    des = load(fullfile(des_center_fp_str.folder, des_center_fp_str.name));
end

end