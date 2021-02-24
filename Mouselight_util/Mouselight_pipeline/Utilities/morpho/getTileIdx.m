function [idxTiles,ia,ic] = getTileIdx(frames,swcData)
%GETTILEIDX Returns id of frames that overlaps with swcData
% 
% [OUTPUTARGS] = GETTILEIDX(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/09/16 14:06:39 $	$Revision: 0.1 $
% Copyright: HHMI 2015

% % create branch objects
% link = double(swcData(:,[1,7]));
% if any(link(:,end)<0)
%     % delete the root node
%     link(link(:,2)<0,:)=[];
% end
% numNodes = max(link(:,1));
%%
% % grab samples from each branch
% if 0
%     L=extractBranches(link,numNodes);
%     subsampledTrace =[];
%     for idx = 1:length(L)
%         set = L(idx).set;
%         if subFactor>0
%             subsampledTrace = [subsampledTrace set(1:round(length(set)/subFactor):end)];
%         else
%             subsampledTrace = [subsampledTrace set(1:end)];
%         end
%     end
%     subsampledTrace = unique(subsampledTrace);
%     subSwcData = swcData(subsampledTrace,:);
% else
%     subSwcData = swcData;
% end
% nSub = length(subsampledTrace);

numNodes = size(swcData,1);
centFrames=squeeze(mean(frames.bbox(:,:,:),1));
numFrames=size(centFrames,2);
L1 = ones(1,numFrames);
hitTable = zeros(1,numNodes);
for idx = 1:numNodes
    locSWC = swcData(idx,3:5);
    dum = (locSWC');
    [minval,minidx] = min(sum((centFrames-dum*L1).^2));
    hitTable(idx) = minidx;
end
[idxTiles,ia,ic] = unique(hitTable,'stable'); % indicies of tiles that are traversed by the swc file
%%



















end
