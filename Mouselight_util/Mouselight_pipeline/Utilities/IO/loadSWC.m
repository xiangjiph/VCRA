function [C,offset,color,header] = loadSWC(swcfile)
%LOADSWC reads a JW formatted swc file and return color/header and data
%fields
%
% [OUTPUTARGS] = LOADSWC(INPUTARGS) Explain usage here
%
% Examples:
%     [swcData,offset,color, header] = loadSWC(swcfile);
%     swcData(:,3:5) = swcData(:,3:5) + ones(size(swcData,1),1)*offset;
%     swcData(:,3:5) = swcData(:,3:5)*scale;
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/14 12:09:52 $	$Revision: 0.1 $
% Copyright: HHMI 2015

%% parse a swc file
[color,offset] = deal([0 0 0]);
header = cell(1);

fid = fopen(swcfile);
tline = fgets(fid);
skipline = 0;
t = 1;
l = 0;
while ischar(tline)
    % assumes that all header info is at the top of the swc file
    if strcmp(tline(1), '#')
        skipline = skipline + 1;
        header{t} = tline;
        t = t+1;
    else
        l=l+1;
        tline = fgets(fid);
        continue;
        %break
    end
    if nargout>1 & strcmp(tline(1:9),'# OFFSET ')
        offset =cellfun(@str2double,strsplit(deblank(tline(10:end))));
    elseif nargout>2 & strcmp(tline(1:8),'# COLOR ')
        color =cellfun(@str2double,strsplit(deblank(tline(9:end)),','));
    end
    tline = fgets(fid);
end
fclose(fid);

fid = fopen(swcfile);
for i=1:skipline
    tline = fgets(fid);
end
C = zeros(l,7);
tline = fgets(fid);
tl = 1;
while ischar(tline)
    C(tl,:) = str2num(tline);
    tl = tl+1;
    tline = fgets(fid);
end
fclose(fid);
end
