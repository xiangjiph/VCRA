function um = pix2um(out,pix)

%%
OX = out.ox;
OY = out.oy;
OZ = out.oz;
OO = [OX OY OZ]/1e3; % in um

SX = out.sx;
SY = out.sy;
SZ = out.sz;
OS = [SX SY SZ]/2^(out.level)/1e3;
%%
pix = pix-1+.5;
%%
um = [pix.*(ones(size(pix,1),1)*OS) + ones(size(pix,1),1)*OO];