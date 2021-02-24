function [out] = configparser(configfile)
%CONFIGPARSER reads a config file and return a structure array based on the
%fields of the configfile
%
% [OUTPUTARGS] = CONFIGPARSER(configfile) 
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% See also: 

% $Author: base $	$Date: 2016/01/14 15:55:07 $	$Revision: 0.1 $
% Copyright: HHMI 2016

if nargin<1
    configfile = 'myparam'
end

fid = fopen(configfile);
out = [];
tline = fgetl(fid);
while ischar(tline)
    if isempty(tline)
    else
        out = checkcase(tline,out);
    end
    tline = fgetl(fid);
end
fclose(fid);

gt=fieldnames(out);
for ii=1:length(gt)
    fieldtxt = strtrim(gt{ii});
    txt = strtrim(out.(gt{ii}));
    if strcmp(gt{ii},'HEADER')
    elseif strcmp(gt{ii},'Tform')
        out.(fieldtxt) = cellfun(@str2num,txt);
    else
        if iscell(txt)
            out.(fieldtxt) = cellfun(@str2num,txt);
        elseif txt(1)=='''' % string, most likely path
            out.(fieldtxt) = eval(txt);
        else
            out.(fieldtxt) = str2num(txt);
        end
    end
end

end

function out = checkcase(tline,out)
if tline(1)=='#' % header skip
    if isfield(out,'HEADER')
        out.HEADER = sprintf('%s%s\n',out.HEADER,tline);
    else
        out.HEADER = sprintf('%s\n',tline);
    end
else
    fd = strsplit(tline,'=');
    if length(fd)<2
        % check for : as splitter
        fd = strsplit(tline,':');
    end
    if length(fd)<2
        error('Couldnt find any argument with "=" or ":"')
    end
    fd_sub = strsplit(strtrim(fd{2}),' ');
    if length(fd_sub)>2
        out.(deblank(fd{1})) = fd_sub;
    else
        out.(deblank(fd{1})) = fd{2};
    end
end

end


