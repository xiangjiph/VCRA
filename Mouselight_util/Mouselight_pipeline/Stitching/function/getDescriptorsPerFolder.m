function [des] = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch,ext_desc)
%GETDESCRIPTORS creates a joint descrtiptor set for multichannel data. 
%
% [OUTPUTARGS] = GETDESCRIPTORS(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
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
    ext_desc = 'txt';
end
myfiles = dir(fullfile(descriptorfolder,scopeloc.relativepaths{10008},['*.',ext_desc]));
fid = fopen(fullfile(myfiles(1).folder,myfiles(1).name));
tLines = fgets(fid);
delimiter = ' '; %or whatever
if length(tLines)>1
    numCols = numel(strfind(tLines,delimiter)) + 1;
    fclose(fid);
else
    error('EMPTY DESCRIPTOR FILE !!')
end
format = repmat('%f ',1,numCols);format(end)=[];
%%
% read descriptor files
numtiles = length(scopeloc.filepath);
desc = [];
parfor_progress(numtiles);
parfor ii=1:numtiles
    parfor_progress;
    for idxch = 1:length(desc_ch)
        %%read descriptors
        %[aa,bb,cc] = fileparts(scopeloc.relativepaths{ii});
        descfile = dir(fullfile(descriptorfolder,scopeloc.relativepaths{ii},['*',desc_ch{idxch},'.',ext_desc]));
        
        if isempty(descfile)
            desc(ii).valid(idxch) = 0;
            desc(ii).value{idxch} = [];
        else
            % load file
            myfid1 = fopen(fullfile(descfile.folder,descfile.name));
            data = textscan(myfid1,format); % round locations, there is a bug in estimation that shifts z-locations 0.5 pix. rounding results in better MSE
            fclose(myfid1);
            if size(data,2)==3
                des0 = [data{1} data{2} data{3}];
            else
                des0 = [[data{1} data{2} data{3}] data{4:numCols}];
            end
            desc(ii).valid(idxch) = 1;
            desc(ii).value{idxch} = des0;
        end
    end
end
parfor_progress(0);

des = cell(1,numtiles);
if length(desc_ch)==1
    for ii=1:numtiles
        des{ii} = desc(ii).value{1};
    end
else
    parfor ii=1:numtiles
        if ~rem(ii,1000)
            ii
        end
        if isempty(desc(ii).value{1})|isempty(desc(ii).value{2})
            continue
        end
        pd2 = pdist2(desc(ii).value{1},desc(ii).value{2});
        [aa1,bb1] = min(pd2,[],1);
        [aa2,bb2] = min(pd2',[],1);% instead of min(pd2,[],2), make it row vector with transpose to prevent dimension error for single row entities
        bb1(aa1>1) = 1;
        bb2(aa2>1) = 1;
        keepthese = [1:length(bb1)]==bb2(bb1);
        des{ii} = desc(ii).value{2}(keepthese,:);
    end
end


end

