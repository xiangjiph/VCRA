function [H,tr] = plot_swc(swcfile,scale)
%PLOT_SWC Wrapper function for TREEs toolbox
%
% [OUTPUTARGS] = PLOT_SWC(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/14 10:45:09 $	$Revision: 0.1 $
% Copyright: HHMI 2015
if nargin <2
    scale = 0;
end
if ischar(swcfile)
    [~,~,EXT] = fileparts(swcfile);
    if strcmp(EXT,'.swc')
        tr=load_tree(swcfile);
    end
elseif isfield(swcfile,'D')
    tr = swcfile;
else
    error('input can be a *.swc file or an intree object')
end
% scale the radius wrto ROI
XYZ = [tr.X tr.Y tr.Z];
% resolution is set to 1/100 of the field
if scale
    rXYZ = min(round(max(XYZ)-min(XYZ)))/200;
else
    rXYZ = 1;
end
tr.D = tr.D*rXYZ;
H = plot_tree(tr,BO_tree(tr));
tr.D = tr.D/rXYZ;
end
