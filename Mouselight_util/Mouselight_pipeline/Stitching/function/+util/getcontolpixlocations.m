function [start_sub,end_sub] = getcontolpixlocations(scopeloc,params,scopeparams)

dims = params.imagesize;
% Q: what's this? 
targetoverlap_um = [5 5 5]; %in um 
N = params.Ndivs;
%% Compute the overlap size in each direction
s1 = median(scopeloc.gridix(:,1:3)); 
% grid ind of the middle tile
i1 = find(scopeloc.gridix(:,1)==s1(1) & scopeloc.gridix(:,2)==s1(2) & scopeloc.gridix(:,3)==s1(3));
s1 = s1+1;
% gird ind of the tiles at the (+1, +1, +1) position w.r.t. middle tile
i2 = find(scopeloc.gridix(:,1)==s1(1) & scopeloc.gridix(:,2)==s1(2) & scopeloc.gridix(:,3)==s1(3));

sdiff = abs(diff(scopeloc.loc([i2 i1],:)))*1000;
overlap_um = round(scopeparams(1).imsize_um - sdiff);
%%
pixsize = scopeparams(1).imsize_um./[scopeparams(1).dims-1];
ovelap_px = round(overlap_um ./ pixsize);
% Number of overlapping pixels in each direction 
targetoverlap_pix = round(targetoverlap_um./pixsize); %in um

start_sub = round(ovelap_px/2 - targetoverlap_pix/2);
% find nearest integer that makes ed-st divisible by number of segments (N)
end_sub = start_sub + floor((dims - 2 * start_sub)/N) * N;

end