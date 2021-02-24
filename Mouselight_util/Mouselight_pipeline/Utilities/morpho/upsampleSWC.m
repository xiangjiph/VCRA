function [upDATA,swcData] = upsampleSWC(swcfile,spacing)
[swcData,offset,color, header] = loadSWC(swcfile);
swcData(:,3:5) = swcData(:,3:5) + ones(size(swcData,1),1)*offset;
swcData(:,3:5) = swcData(:,3:5)*1000; % in nm

%%
upDATA = [];
% upsample points that are consequtive
upthese = find(swcData(:,7)~=-1);
updata_ = [];
for jj=upthese(:)'
    st = swcData(jj,3:5);
    en = swcData(swcData(jj,7),3:5);
    sl = en-st;
    pd = norm(sl);
    sp = [[1:spacing:pd-1 pd]/pd];
    %updata = [updata sl(:)*sp+st(:)*ones(1,length(sp))];
    updata_{jj} = [sl(:)*sp+st(:)*ones(1,length(sp))];
end
upDATA = cat(1,upDATA,[updata_{:}]');

