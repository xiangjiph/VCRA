function [brainSize,RR,chunk_dims,rank] = h5parser(myh5,myh5prob)
%H5PARSER retrive related fields from a large h5 file. using h5info is
%exteremly slow, especialy if the file is hosted on a shared drive
%
% [OUTPUTARGS] = H5PARSER(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/08/03 12:47:24 $	$Revision: 0.1 $
% Copyright: HHMI 2017

fid = H5F.open(myh5);
dset_id = H5D.open(fid,myh5prob);
space = H5D.get_space(dset_id);
[~,dims] = H5S.get_simple_extent_dims(space);
H5S.close(space);

dcpl = H5D.get_create_plist(dset_id);
[rank,chunk_dims] = H5P.get_chunk(dcpl);
H5P.close(dcpl);

brainSize = dims([3 2 1]);
chunk_dims = chunk_dims([3 2 1]);
RR = h5read(myh5,[myh5prob,'_props/ROI']);

H5D.close(dset_id);
H5F.close(fid);


end
