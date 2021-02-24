function [projLoc] = patch2render(labs,idx,frames,scale)
%PATCH2RENDER Maps tile id back to rendered location
% 
% [OUTPUTARGS] = PATCH2RENDER(INPUTARGS) Explain usage here
% 
% Inputs: 
%   labs: tile label matrix
%   idx: tile id
%   frames: yaml structure
%   scale: unit conversion
% 
% Outputs: 
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/08/31 11:52:31 $	$Revision: 0.1 $
% Copyright: HHMI 2015

sepidx = findstr(labs{idx},labs{idx}(1));
tilename1 = labs{idx}(sepidx(1)+1:sepidx(2)-1);
% find subfolderidx
val1 = ~cellfun(@isempty,strfind(frames.path,tilename1));

strs = strsplit(labs{idx},' ');
patchloc = [str2double(strs{end-1}) str2double(strs{end})];
tilename2 = labs{idx}(sepidx(2)+1:sepidx(2)+5);
val2 = ~cellfun(@isempty,strfind(frames.path,strrep(strs{1},'\','/')));
tileidx = find(val1&val2);
% disp(['Found:', frames.path{tileidx}, 'at idx:',num2str(tileidx)]);
Tform = frames.tform(:,:,tileidx);
projLoc = [Tform(1:3,:)*[patchloc 0 1]'/scale]';



end
