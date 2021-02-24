function desc = truncateDesc(desc,maxnumofdesc,thr)
if nargin<2
    maxnumofdesc = 10e3;
    thr = 0.1;
elseif nargin<3
    thr = 0.1;
end

if maxnumofdesc<size(desc,1)
    [vals,indssorted]=sort(desc(:,5),'descend'); % sort based on raw intensity
    validinds = indssorted(1:min(maxnumofdesc,length(indssorted)));
    desc = double(desc(validinds,:));
else
    desc = double(desc(desc(:,5)>thr,:));
    desc(desc(:,3)>249,:)=[];
end
