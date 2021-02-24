function [maxval,maxarg,dists] = dist2line(x1,x2,x0)
%DIST2LINE finds the distance of x0 to line specified by x1&x2
% 
% [OUTPUTARGS] = DIST2LINE(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/08/19 10:30:44 $	$Revision: 0.1 $
% Copyright: HHMI 2015
% x2-x1)*(y1-y0)-(x1-x0)*(y2-y1)//sqrt((x2-x1)^2+
% d = abs((x2(1,:)-x1(1,:)) * (x1(2)-x0(2)) - (x1(1,:)-x0(1,:)).*(x2(2,:)-x1(2,:)))/norm(x2-x1);

m10 = (x1(2)-x0(2,:))./(x1(1)-x0(1,:));
m12 = (x1(2)-x2(2))/(x1(1)-x2(1));
dists = abs((m10-m12).*(x2(1)-x1(1)).*(x1(1)-x0(1,:))/norm(x2-x1));


[maxval,maxloc] = max(dists);
maxarg = x0(:,maxloc);
end
