function [outtree] = smoothbranches(in_)
%SMOOTHBRANCHES Gets an affinity matrix and node locations and return a
%smoothed locations for each segments
% 
% [OUTPUTARGS] = SMOOTHBRANCHES(INPUTARGS) Explain usage here
% 
% Inputs: 
% 
% Outputs: 
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/10/21 16:33:23 $	$Revision: 0.1 $
% Copyright: HHMI 2015

%% smooth
clc
A = in_.dA;
XYZ = [in_.X in_.Y in_.Z];
R = in_.R;
D = in_.D;
%%
% ani = [1 1 3];
% var = sum(ani.^2);

thr_branchlen = 3;
[L,list] = getBranches(A);
XYZori = XYZ;
Lori = L;

for ii=1:length(L)
    %%
    set_ii = L(ii).set;
    if length(set_ii)>thr_branchlen % ong enough branch
        % smoothing with running avarage
        xyz = XYZ(set_ii,:);
        for jj=3 % only on z
            X = xyz(:,jj);
            XYZ(set_ii,jj) = medfilt1(medfilt1(X,thr_branchlen),thr_branchlen);
        end
        % down sampling
        [out] = downSamplePath(XYZ(set_ii,:),5); 
        % update branch
        L(ii).set =set_ii(out);
    end
end

%%
% build a new affinity graph 
% [A_,XYZ_] = branch2tree(L,A,XYZout);
[A_,idx] = branch2tree(L);
% [L2,list2] = getBranches(A_);

XYZ_ = XYZ(idx,:);
R_ = R(idx,:);
D_ = D(idx,:);
outtree.dA = A_;
outtree.X = XYZ_(:,1);
outtree.Y = XYZ_(:,2);
outtree.Z = XYZ_(:,3);
outtree.R = R_;
outtree.D = D_;
% %%
% 
% figure_(100);
% cla
% hold on
% gplot3(A,XYZ,'--b','LineWidth',2);
% % gplot3(A,XYZout,'r','LineWidth',3);
% gplot3(A_,XYZ_,'g','LineWidth',3);
% % axis equal tight
% 
% %%
% % close all
% % figure_,
% cla
% hold on
% gplot3(A_,XYZ_,'r','LineWidth',2)
% gplot3(A,XYZ,'--b')
% axis equal tight
% %%
% figure_,
% hold on
% gplot3(A,XYZ,'r','LineWidth',2)
% gplot3(A(set_ii,set_ii),XYZ(set_ii,:),'b','LineWidth',2)
% axis equal tight

end
