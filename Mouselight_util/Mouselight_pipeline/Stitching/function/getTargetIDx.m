function targetidx = getTargetIDx(scopeloc,neighbors)
if 1
    %%
    if 1
        st=[192,128,2133]-1;%[192,124,1049];
        ed=st+2;
    elseif 0
        for ii=1:1e6, if length(strfind(scopeloc.filepath{ii},'/2017-04-26/01/01120')),break,end,end,ii
        st=scopeloc.gridix(ii,1:3)-1;
        ed=st+[1 1 1]*2;
    elseif 0
        % oct12
        st = [1 1 1]
        ed = [35 9 40]
        
    end
    targetidx = find(scopeloc.gridix(:,1)>=st(1)&scopeloc.gridix(:,1)<=ed(1)&...
        scopeloc.gridix(:,2)>=st(2)&scopeloc.gridix(:,2)<=ed(2)&...
        scopeloc.gridix(:,3)>=st(3)&scopeloc.gridix(:,3)<=ed(3))
%     targetidx = [1971,neighbors(1971,5)]
    %         targetidx=[6962 6991]
elseif 0
    inds_ = 568;%inds(1)';
    neigs = neighbors(inds_,checkthese);
    targetidx = neigs(1:2)
else
    if 1
        med = scopeloc.gridix(9512,1:3)
    else
        med = median(scopeloc.gridix(scopeloc.gridix(:,3)==120,:))
    end
    targetidx = [];
    % med = scopeloc.gridix(12300,:)
    params.outfile = sprintf('%sxyFCROI-%d%d%d_%s.control.yml',experimentfolder,med,date);
    for t = med(3)-1:med(3)+1%latticeZRange(40)'%22:24%
        %%
        disp(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
        ix = (scopeloc.gridix(:,3)==t)&(scopeloc.gridix(:,2)>=(med(2)-1)&scopeloc.gridix(:,2)<=(med(2)+1))&...
            (scopeloc.gridix(:,1)>=(med(1)-1)&scopeloc.gridix(:,1)<=(med(1)+1));
        targetidx = [targetidx find(ix)]
    end
end