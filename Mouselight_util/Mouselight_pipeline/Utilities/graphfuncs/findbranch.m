function br=findbranch(A,X)
algorithmkeys = {'spb','dij','bel','spa'};
algorithm = 2;
debug_level = 0;
directed = 0;
smoothpath=1;
sizethreshold = 5;

Ain=A;
leafnodes = find(sum(A,2)==1);
regnodes = find(sum(A,2)==2);
juncnodes = find(sum(A,2)>2);
cricset=union(leafnodes,juncnodes);
% delete junction node edges
A(juncnodes,:) = 0;
A(:,juncnodes) = 0;

CompsC = conncomp(graph(A),'OutputForm','cell');
S = length(CompsC);
Y = cellfun(@length,CompsC);
%%
tipids = zeros(S,2);
br=[];
iter=1;
for ii=1:S
    if S>10 & ~rem(ii,round(S/10))
%         ii
    end
    inds = CompsC{ii};
    if length(inds())<2
        % skip
        %         break
        continue
    end
    A_ = A(inds,inds);
    subs_ = X(inds,:);
    tips = find(sum(A_,1)==1);
    %%
    % for each tip find neighbor
    if 1
        indtips = inds(tips);
        [aa,bb]=ind2sub([size(Ain,1) length(indtips)],find(Ain(:,indtips)));
        % [aa indtips(bb)']
        indnodes = union(aa,inds);
        %%%%%
        %tar1 = find(Ain(:,indtips(1)));
        % check if they are junction/leaf
    else
        indnodes = inds;
    end
%     if length(indnodes)==3 % > leaf or self loop at tip
%         
%         
%         
%     end
    %%
    A_ = Ain(indnodes,indnodes);
    subs_ = X(indnodes,:);
    extendedtips = find(sum(A_,1)==1);
    if isempty(extendedtips) % small tip loop
        % skip tip loops
        continue
    end
    if version('-release')=='2015b'
        [dist,pred] = graphalgs(algorithmkeys{algorithm},...
            debug_level,directed,A_,extendedtips(1));
        
        path = graphpred2path(pred,extendedtips(2));
    else
        G_=graph(A_);path = shortestpath(G_,extendedtips(1),extendedtips(2));
    end
    
    set = indnodes(path);
    X_ = subs_(path,:);
    if smoothpath & length(indnodes)>sizethreshold
        for jj=1:3
            X_(2:end-1,jj) = medfilt1(medfilt1(X_(2:end-1,jj),sizethreshold),sizethreshold);
        end
    end
    br(iter).termnodes = indnodes(extendedtips);
    br(iter).set = indnodes(path);
    br(iter).subs = X_;
    iter=iter+1;
end
%%
