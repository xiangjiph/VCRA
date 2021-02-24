function D = calcDists(branches,connG)

nA = length(branches);
[ix,iy] = find(connG);
num_edges = length(ix);

%%
% clc
warning off
w=zeros(num_edges,4);
conn=zeros(num_edges,1);
tic
parfor idx=1:num_edges
    br1 = branches(ix(idx)).inds;%path;
    sub1 = branches(ix(idx)).subs;%subs(br1,:);
    br2 = branches(iy(idx)).inds;%path;
    sub2 = branches(iy(idx)).subs;%subs(br2,:);
    %
    if length(br2)<4|length(br1)<4
        w(idx,:) = [0,10,0,0];
    else
        % align tips
        pd=pdist2(sub1([1 end],:),sub2([1 end],:));
        [aa,bb]=min(pd(:));
        [ax,ay] = ind2sub(size(pd),bb);
        % alignment score
        type = 0;
        if (ax==1&ay==1)
            type = 1;
            % flip one of them
            ds=alignmentS(flipud(br1(:)),flipud(sub1),br2(:),sub2);
        elseif (ax==2&ay==1)
            %
            type = 2;
            ds=alignmentS(br1(:),sub1,br2(:),sub2);
        elseif (ax==1&ay==2)
            %
            type = 3;
            ds=alignmentS(br2(:),sub2,br1(:),sub1);
        else
            type = 4;
            ds=alignmentS(br1(:),sub1,flipud(br2(:)),flipud(sub2));
        end
        w(idx,:) = ds;
        conn(idx,:) = type;
    end
end
toc
warning on

%%
w(w==0)=eps;
clear D
for ii=1:size(w,2)
    D{ii} = sparse(ix,iy,w(:,ii),nA,nA);
end
D{end+1} = sparse(ix,iy,conn(:),nA,nA);

% D(1): Euclidean
% D(2): theta
% D(3): PCA
% D(4): KL
% D(5 [end]): hit

% D{1} = sparse(ix,iy,w(:,1),nA,nA);
% D{2} = sparse(ix,iy,w(:,2),nA,nA);
% D{3} = sparse(ix,iy,w(:,3),nA,nA);
%
%
% W1 = inf(nA);
% W1(connG) = w(:,1); % theta
% W2 = inf(nA);
% W2(connG) = w(:,2); % PCA
% W3 = inf(nA);
% W3(connG) = w(:,3); % KL
% W = cat(3,max(W1,W1'),max(W2,W2'),max(W3,W3'));
end

%%%%%%%
function alignmentdist = alignmentS(br1,sub1,br2,sub2)

%%
alignmentdist = [];
alignmentdist(1) = norm(sub1(end,:)-sub2(1,:));
l1 = length(br1);
l2 = length(br2);

st1 = max(1,l1-9);
ed1 = l1;
st2 = 1;
ed2 = min(10,l2);
kappa=1e-5;
ndims = size(sub1,2);
sub1_ = sub1(st1:ed1,:)+kappa*(rand(ed1-st1+1,ndims)-.5);
sub1_ = (sub1_-ones(size(sub1_,1),1)*mean(sub1_,1));
sub2_ = sub2(st2:ed2,:)+kappa*(rand(ed2-st2+1,ndims)-.5);
sub2_ = (sub2_-ones(size(sub2_,1),1)*mean(sub2_,1));
A = cov(sub1_);
B = cov(sub2_);

[eA,uA] = eig(A);
[eB,uB] = eig(B);
e1 = eA(:,ndims);
e2 = eB(:,ndims);
thet = acos(max(-1,min(1,sum(e1.*e2))));
% alignmentdist(1) = thet;%min(thet,pi-thet);
alignmentdist(end+1) = min(thet,pi-thet);

% projection distance
if 0
    e1 = normr(sub1(end,:)-sub1(1,:))';
    e2 = normr(sub2(end,:)-sub1(2,:))';
    panchor = sub1(end,:);
    v = panchor-mean(sub2(st2:ed2,:),1);
    bet1 = acos(max(-1,min(1,sum(e2.*normc(v(:))))));
    panchor = sub2(1,:);
    v = panchor-mean(sub1(st1:ed1,:),1);
    bet2 = acos(max(-1,min(1,sum(e1.*normc(v(:))))));
    score(3) = sin(max(bet1,bet2));
else
    X = [sub1(st1:ed1,:);sub2(st2:ed2,:)];
    [n,p] = size(X);
    meanX = mean(X,1);
    if 1
        X_ = [X-meanX];
        [coeff,roots] = eig(X_'*X_);
        %[coeff,roots] = eig(cov(X_));
        [roots,inds] = sort(diag(roots),'descend');
        coeff = coeff(:,inds);
        reconstructed_signal = (X_*coeff(:,1))*coeff(:,1)';
        rmse = sqrt(mean(sum((X_-reconstructed_signal).^2,2)));
    else
        [coeff,score,roots] = pca(X);
        Xfit = repmat(meanX,n,1) + score(:,1)*coeff(:,1)';
        rmse = sqrt(mean(sum((X-Xfit).^2,2)));
    end
    alignmentdist(end+1) = rmse;
end
warning('off')
score1 = min(1e3,1/2*(trace(A\B) -length(sub1(1,:)) - log(det(A)/det(B))));
score2 = min(1e5,1/2*(trace(B\A) -length(sub1(1,:)) - log(det(B)/det(A))));
warning('on')
alignmentdist(end+1) = sqrt(abs(score1*score2));
%%
% it=0
% figure(123),cla
% hold on
% % gplot3(A,subs,'b')
% myplot3(sub1(end-[10:-1:0],:),'-r')
% myplot3(sub2(1:10,:),'-m')
end



















