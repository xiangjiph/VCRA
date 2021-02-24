function [nooptim,control_t_bot12,control_tp1_top12,zlim_cntrl] = ...
    nomatchoptim(idxt,params,tileneighbors,afftile,pixstats,zlim_cntrl,...
    corrctrlpnttmp,Fxt,Fyt,Fzt,Fxtp1,Fytp1,Fztp1,zlimdefaults)

nooptim = false;
htop = params.htop;
Nlayer = params.Nlayer;
Npts = (params.Ndivs+1).^2;
dims = params.imagesize;
order = params.order;
if nargin<14
    zlimdefaults = params.zlimdefaults;
end

noopt = 0; % no optimization

idxtm1 = tileneighbors(idxt,6);
idxtp1 = tileneighbors(idxt,7);
tformt = afftile(:,:,idxt);
tstats = pixstats(idxt,:); % [idxt  idxtp1 med[idxt,idxtp1] min[idxt,idxtp1] max[idxt,idxtp1]];
if ~isnan(idxtm1)
    tformtm1 = afftile(:,:,idxtm1);
    tstats_tm1 = pixstats(idxtm1,:); % [idxt  idxtp1 med[idxt,idxtp1] min[idxt,idxtp1] max[idxt,idxtp1]];
end
if ~isnan(idxtp1)
    tformtp1 = afftile(:,:,idxtp1);
    tstats_tp1 = pixstats(idxtp1,:); % [idxt  idxtp1 med[idxt,idxtp1] min[idxt,idxtp1] max[idxt,idxtp1]];
end


% same as above for the bottom of the current tile t.
zlim_3 = zlimdefaults(3);%dims(3)-11;
zlim_4 = zlimdefaults(4);%dims(3)-1;
tstats(8) = max(pixstats(pixstats(:,8)<150,8)); % added a hard cap at 150 to prevent artifacts due to curation overwrite mistakes
tstats(6) = min(pixstats(:,6));
%% Exactly the same as optimpertile line 88 - 106
corrctrlpnttmp(:,3) = zlim_3;
contr_t_bot1 = [corrctrlpnttmp, ones(Npts,1)] * tformt';
vect_bot1 = [Fxt(contr_t_bot1), Fyt(contr_t_bot1), Fzt(contr_t_bot1)];
if noopt
    vect_bot1(:) = 0;
end
contr_t_bot1_shifted = contr_t_bot1 + vect_bot1;

corrctrlpnttmp(:,3) = zlim_4;
contr_t_bot2 = [corrctrlpnttmp ones(Npts,1)] * tformt';
vect_bot2 = [Fxt(contr_t_bot2), Fyt(contr_t_bot2), Fzt(contr_t_bot2)];
if noopt
    vect_bot2(:) = 0;
end
contr_t_bot2_shifted = contr_t_bot2 + vect_bot2;

% book it
control_t_bot12 = [contr_t_bot1_shifted; contr_t_bot2_shifted];
%% %%%%%%%%%%%%%%%%%%%%%
% SEARCH
for zlim_2 = tstats(8)-1:-1:max(1,min(tstats(8)-1,tstats(6))) %
    %zlim_2 = max(zlim_1_bound+1,tstats(8)-hbuf);
    corrctrlpnttmp(:,3)= zlim_2;
    
    contr_tp1_top2 = [corrctrlpnttmp ones(Npts,1)]*tformtp1';
    % shift it
    vectp1_top2=[Fxtp1(contr_tp1_top2) Fytp1(contr_tp1_top2) Fztp1(contr_tp1_top2)];
    if noopt
        vectp1_top2(:) = 0;
    end
    contr_tp1_top2_shifted = contr_tp1_top2+vectp1_top2;
    
    % make sure that there is overlap for all control points
    if all(contr_t_bot2_shifted(:,3)-contr_tp1_top2_shifted(:,3) > 0)
        %zlim_2
        break
    end
end
%%%%%%
%%%%%%
% top1: then get top of overlap on layer t at crop "h"
% (1)
zlim_1 = max(0,zlim_2-htop);
corrctrlpnttmp(:,3)= zlim_1;
contr_tp1_top1 = [corrctrlpnttmp ones(Npts,1)]*tformtp1';
% shift it
vectp1_top1=[Fxtp1(contr_tp1_top1) Fytp1(contr_tp1_top1) Fztp1(contr_tp1_top1)];
if noopt
    vectp1_top1(:) = 0;
end
contr_tp1_top1_shifted = contr_tp1_top1+vectp1_top1;
control_tp1_top12 = [contr_tp1_top1_shifted;contr_tp1_top2_shifted];
zlim_cntrl(1:2,idxtp1) = [zlim_1;zlim_2];
zlim_cntrl(3:4,idxt) = [zlim_3;zlim_4];
% check if zlim ordering changes due to knn search
if zlim_cntrl(2,idxt) >= zlim_cntrl(3,idxt)
    if zlim_cntrl(1,idxt) < zlim_cntrl(3,idxt)-2 % make sure there is enough space with the first control point
        zlim_cntrl(2,idxt) = zlim_cntrl(3,idxt)-1;
    else
        if zlim_cntrl(3,idxt)>3
            zlim_cntrl(2,idxt) = zlim_cntrl(3,idxt)-1;
            zlim_cntrl(1,idxt) = zlim_cntrl(3,idxt)-2;
        end
    end
end


