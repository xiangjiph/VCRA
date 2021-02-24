function [paireddescriptor, scopeparams, R, curvemodel,scopeparams_, paireddescriptor_,curvemodel_] = ...
    estimateFCpertile(descriptors,neighbors,scopeloc,params)
%ESTIMATEFCPERTILE Summary of this function goes here
% 
% [OUTPUTARGS] = ESTIMATEFCPERTILE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/25 10:18:10 $	$Revision: 0.1 $
% Copyright: HHMI 2017
%%
checkthese = [1 4 5 7]; % 0 - right - bottom - below
imsize_um = params.imsize_um;
mkdir('./matfiles')
if 1
    [paireddescriptor,R,curvemodel] = xymatch_old(...
        descriptors,neighbors(:,checkthese),scopeloc,params);
    save ./matfiles/xypaireddescriptor paireddescriptor R curvemodel
else
    load ./matfiles/xypaireddescriptor paireddescriptor R curvemodel
end
%%
if 1
    if 1
        % based on 4 neig + z regularizor
        [scopeparams,scopeparams_,paireddescriptor_,curvemodel_] = homographyPerTile6Neighbor(...
            params,neighbors,scopeloc,paireddescriptor,R,curvemodel);
        % based on 2 neig + z regularizor
%         [scopeparams,scopeparams_,paireddescriptor_,curvemodel_] = homographyPerTile(...
%             beadparams,neighbors(:,checkthese),scopeloc,paireddescriptor,R,curvemodel,imsize_um);
        save ./matfiles/scopeparams scopeparams scopeparams_ paireddescriptor_ curvemodel_
    else
        load ./matfiles/scopeparams scopeparams scopeparams_ paireddescriptor_ curvemodel_
    end
end
end
