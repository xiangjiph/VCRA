function nummatches = matchLength(featmap,directions)
numtiles = length(featmap);
nummatches = zeros(1,length(featmap));
for ii=1:numtiles
    feat = featmap(ii).(genvarname(directions));
    if isempty(feat)
    else
        nummatches(ii) = size(feat.paireddescriptor.X,1);
    end
end