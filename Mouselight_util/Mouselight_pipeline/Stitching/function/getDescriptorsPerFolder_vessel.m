function [des] = getDescriptorsPerFolder_vessel(descriptorfolder,scopeloc,desc_ch,ext_desc)
%GETDESCRIPTORS creates a joint descrtiptor set for multichannel data. 
%
% [OUTPUTARGS] = GETDESCRIPTORS(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%   des: N-by-1 cell array, N is the number of tiles in the folder 
%       each cell contains the descriptor structure
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/21 16:41:34 $	$Revision: 0.1 $
% Copyright: HHMI 2016

%% get a sample file to check format
if nargin<4
    ext_desc = 'mat';
end
num_tiles = length(scopeloc.filepath);
des = cell(1,num_tiles);
% If we have multiple channel descriptors in terms of MATLAB stucture,
% merge the fields directly. Notice that 
parfor_progress(num_tiles);
parfor iter_tile = 1 : num_tiles
    parfor_progress;
    tmp_merged_struct = struct;
    for idxch = 1:length(desc_ch)
        descfile = dir(fullfile(descriptorfolder,scopeloc.relativepaths{iter_tile},['*desc.', desc_ch{idxch},'.',ext_desc]));
        tmp_struct = load(fullfile(descfile.folder, descfile.name));
        tmp_field_name = fieldnames(tmp_struct);
        for iter_field = 1 : numel(tmp_field_name)
            if ~isfield(tmp_merged_struct, strip(tmp_field_name{iter_field}))
                new_field_name = strip(tmp_field_name{iter_field});
            else
                new_field_name = sprintf('Ch_%s_%s', desc_ch{idxch}, strip(tmp_field_name{iter_field}));
                warning('The filed name alreadly exist. Rename the field according to channel');
            end
            tmp_merged_struct.(new_field_name) = tmp_struct.(tmp_field_name{iter_field});
        end
    end
    des{iter_tile} = tmp_merged_struct;
end
parfor_progress(0);
%% For feature merging and selection 
% % read descriptor files
% desc = [];
% parfor_progress(num_tiles);
% parfor ii=1:num_tiles
%     parfor_progress;
%     for idxch = 1:length(desc_ch)
%         %%read descriptors
%         %[aa,bb,cc] = fileparts(scopeloc.relativepaths{ii});
%         descfile = dir(fullfile(descriptorfolder,scopeloc.relativepaths{ii},['*',desc_ch{idxch},'.',ext_desc]));
%         
%         if isempty(descfile)
%             desc(ii).valid(idxch) = 0;
%             desc(ii).value{idxch} = [];
%         else
%             % load file
%             myfid1 = fopen(fullfile(descfile.folder,descfile.name));
%             data = textscan(myfid1,format); % round locations, there is a bug in estimation that shifts z-locations 0.5 pix. rounding results in better MSE
%             fclose(myfid1);
%             if size(data,2)==3
%                 des0 = [data{1} data{2} data{3}];
%             else
%                 des0 = [[data{1} data{2} data{3}] data{4:numCols}];
%             end
%             desc(ii).valid(idxch) = 1;
%             desc(ii).value{idxch} = des0;
%         end
%     end
% end
% parfor_progress(0);
% 
% if length(desc_ch)==1
%     for ii=1:num_tiles
%         des{ii} = desc(ii).value{1};
%     end
% else
%     parfor ii=1:num_tiles
%         if ~rem(ii,1000)
%             ii
%         end
%         if isempty(desc(ii).value{1})|isempty(desc(ii).value{2})
%             continue
%         end
%         pd2 = pdist2(desc(ii).value{1},desc(ii).value{2});
%         [aa1,bb1] = min(pd2,[],1);
%         [aa2,bb2] = min(pd2',[],1);% instead of min(pd2,[],2), make it row vector with transpose to prevent dimension error for single row entities
%         bb1(aa1>1) = 1;
%         bb2(aa2>1) = 1;
%         keepthese = [1:length(bb1)]==bb2(bb1);
%         des{ii} = desc(ii).value{2}(keepthese,:);
%     end
% end


end

