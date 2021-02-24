function [threshold,x1,x2,maxarg,vals,bins] = get1DThresh(In,nbins,perc)
%GETTHRESH finds the binary threshold based on histogram maximal distance
% 
% [OUTPUTARGS] = GETTHRESH(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/08/19 10:37:18 $	$Revision: 0.1 $
% Copyright: HHMI 2015
%%
if nargin<2
    nbins = 256;
    perc = 0;
elseif nargin<3
    perc=0;
end
[vals,bins] = hist(double(In(:)),nbins);

[x12,locmax] = max(vals);

% below works better if there are multiple peaks (dark pixel and background have different peaks)
% append 0 to both sides
xleft = [0 vals(1:end-1)];
xright = [vals(2:end) 0];
maximas = find(vals(1:nbins/2)>xleft(1:nbins/2) &vals(1:nbins/2)>xright(1:nbins/2) & vals(1:nbins/2)>vals(locmax)/2);
[peakLoc] = util.peakfinder(vals);
locmax = peakLoc(1);
x12 = vals(locmax);

x1 = [bins(locmax) x12]';

% % apply suppresion
% [vals,bins] = hist(double(In(In>bins(locmax))),100);
if perc
% find the percentile that has %95 of right hand data
    idx2 = find((cumsum(vals)/sum(vals(:)))>perc,1,'first');
else
    idx2 = length(bins);
end
x2 = [bins(idx2) 0]';
x0 = [bins(locmax+1:end);vals(locmax+1:end)] ;
% make sure solution is in the convex region (this is necessary for perc 
% calculations)
x0 = x0(:,x0(1,:)<x2(1) & x0(1,:)>x1(1));
[maxval,maxarg,d]=dist2line(x1,x2,x0);
if numel(maxarg)
    threshold = maxarg(1);
else
    maxIn = max(In(:));
    threshold = max(1,graythresh(In/maxIn)*maxIn); % for heavy peaked distributions, OTSU returns 0
end
end
