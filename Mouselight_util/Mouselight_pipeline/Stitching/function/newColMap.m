function newMap = newColMap(centerPoint,scalingIntensity,dataMin,dataMax);
if nargin<2
    newMap = scalemap(centerPoint);
    return
end
if nargin<3
    centerPoint = 25;
    scalingIntensity = 5;
    dataMin = 0;
    dataMax = 255;
end
numbin = 2^9;
cMap = redbluecmap(256);
x = 1:length(cMap);
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*(numbin-1)/max(x)+1;
newMap = interp1(x, cMap, 1:numbin);
newMap = scalemap(numbin);
end

function newMap = scalemap(m)
cMap = redbluecmap(12);
newMap = interp1(linspace(-1,1,11), cMap, linspace(-1,1,m));
end