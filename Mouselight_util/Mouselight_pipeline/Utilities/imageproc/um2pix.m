function pix = um2pix(params,um)
% converts um coordinates onto voxels locations in hdf5 file
OX = params.ox;
OY = params.oy;
OZ = params.oz;
OO = [OX OY OZ]/1e3; % in um

SX = params.sx;
SY = params.sy;
SZ = params.sz;
level = params.level;
OS = [SX SY SZ]/2^(level)/1e3;
%%
pix = [um - ones(size(um,1),1)*OO]./(ones(size(um,1),1)*OS);
pix = round(pix + .5);

