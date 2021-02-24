function mylist = parfor_solver(Idx,dist,nodeBrid)
% splits large data to be processed in parfor
numbr = length(Idx)/2;
mylist = cell(numbr,1);
neig_size=cellfun(@length,nodeBrid);
%%
try parfor_progress(0);catch;end
parfor_progress(numbr)
tic
for ibr=1:numbr
    %%
    from = ibr;
    inds=[Idx{(ibr-1)*2+1:ibr*2}]; % all nodes around two tips
    dists=[dist{(ibr-1)*2+1:ibr*2}]; % distance to querry branch of a node
    a_ = cell(1,length(inds));
    for jind = 1:length(inds)
        brid = nodeBrid{inds(jind)};
        if isempty(brid)
        elseif length(brid)>1
            a_{jind} = [ones(length(brid),1)*from brid(:) ones(length(brid),1)*dists(jind)];
        else
            a_{jind} = [from brid dists(jind)];
        end
    end
    a_ = cat(1,a_{:});
    [dx,ix] = sort(a_(:,3));
    [to,IA,IC] = unique(a_(ix,2)); % branch ids of all nodes
    % [to dx(IA) wq(:)]
    dfromto = dx(IA);
    dfromto(to==from)=[];
    dfromto(dfromto==0) = eps;
    to(to==from)=[];
    mylist{ibr} = [from*ones(length(to),1) to(:) dfromto(:)];
    parfor_progress;
end
parfor_progress(0)
sprintf('Branch-Con DONE in %d sec',round(toc))



