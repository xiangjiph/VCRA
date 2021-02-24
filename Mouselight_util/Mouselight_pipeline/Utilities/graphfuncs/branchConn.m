function distBr = branchConn(branches,subs,nodeBrid,querdist)

numbr = length(branches);
querrylocs = zeros(numbr*2,3);
for ibr = 1:numbr
    querrylocs((ibr-1)*2+1:ibr*2,:) = branches(ibr).subs([1 end],:);
end
%%
% build kdtree around branch tip points
[Idx,dist] = rangesearch(subs,querrylocs,querdist);
%%
mylist = graphfuncs.parfor_solver(Idx,dist,nodeBrid);
%%
mylist = cat(1,mylist{:});
%%
% create branch list/lookup
distBr = sparse(mylist(:,1),mylist(:,2),mylist(:,3),numbr,numbr);
