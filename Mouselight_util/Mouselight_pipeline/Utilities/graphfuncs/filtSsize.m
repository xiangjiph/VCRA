function [G,subs] = filtSsize(G,subs,argin)
%FILTSSIZE Summary of this function goes here
% 
% [OUTPUTARGS] = FILTGSIZE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/05 16:35:55 $	$Revision: 0.1 $
% Copyright: HHMI 2017

sizeThr = argin{1};
opt = argin{2};
% get h5 files
h5_file = opt.inputh5;
h5_ch1 = strrep(h5_file,'prob0','prob1');
h5_ch0 = strrep(h5_ch1,'prob1','prob0');

CompsC = conncomp(G,'OutputForm','cell');
% S = length(CompsC);
Y = cellfun(@length,CompsC);
deletethese = (Y<sizeThr);
cpl = [CompsC{deletethese}];
%%
delsubs = subs(cpl,:);
%%
%parcel
%%
af = zeros(size(delsubs));
bb=[3 3 3]
grapvol1 = zeros([bb 1e3]);
grapvol2 = zeros([bb 1e3]);
tic
for ii=1:1e3%size(delsubs,1)
    starts = delsubs(ii,:)-1;
    grapvol1(:,:,:,ii) = double(squeeze(h5read(h5_ch0,'/prob0',starts,bb)));
end
toc
tic
for ii=1:1e3%size(delsubs,1)
    starts = delsubs(ii,:)-1;
    grapvol2(:,:,:,ii) = double(squeeze(h5read(h5_ch1,'/prob1',starts,bb)));
end
toc
grapvol = grapvol1.*grapvol2;
%%
    grapvol2 = grapvol1;
%     grapvol2 = double(squeeze(h5read(h5_ch1,'/prob1',starts,[5 5 5])));
    gv = grapvol1.*grapvol2;
    af(ii) = sum(gv(:))/norm(grapvol1(:))/norm(grapvol2(:))>.5;
    

%%
% check if they exist in both channel
for ii=1:deletethese
    CompsC{deletethese(ii)}
    
end
%%
%%
G = rmnode(G,cpl);
subs(cpl,:)=[];
end
