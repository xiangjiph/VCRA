function [indicies] = sampleXum(cdist,sampling,set_ii)
%SAMPLEXUM Summary of this function goes here
% 
% [OUTPUTARGS] = SAMPLEXUM(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/12/22 16:42:50 $	$Revision: 0.1 $
% Copyright: HHMI 2015

filterrad = 5; % fragments are trimmed around each branch for cleaned up visualization
start = set_ii(find(cdist>min(sampling,filterrad),1)); % find the first index that is 5um away
indicies = zeros(1,ceil(cdist(end)/sampling));
for jj=1:ceil(cdist(end)/sampling)
    idx = find(cdist>jj*sampling,1);
    if isempty(idx)
    else
        indicies(jj) = set_ii(idx);
    end
end
if ~(indicies(end))
    indicies(end) = set_ii(end);
end
indicies = unique([start indicies]','stable');

end
