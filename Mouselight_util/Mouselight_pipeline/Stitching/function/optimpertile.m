function [nooptim,control_t_bot12,control_tp1_top12,zlim_cntrl]=...
    optimpertile(idxt,params,tileneighbors,afftile,pixstats,zlim_cntrl,corr_ctrl_pnt_tmp,Fxt,Fyt,Fzt,Fxtp1,Fytp1,Fztp1)
%% Parameters
nooptim = false;
[control_t_bot12,control_tp1_top12] = deal([]);
htop = params.htop;
Nlayer = params.Nlayer;
Npts = (params.Ndivs+1).^2;
dims = params.imagesize;
order = params.order;
zlimdefaults = params.zlimdefaults;

noopt = 0; % no optimization

idx_tm1 = tileneighbors(idxt,6); % sub_3 - 1
idx_tp1 = tileneighbors(idxt,7); % sub_3 + 1
tformt = afftile(:,:,idxt);
tstats = pixstats(idxt,:); % [idxt  idxtp1 med[idxt,idxtp1] min[idxt,idxtp1] max[idxt,idxtp1]];
if ~isnan(idx_tm1)
    tformtm1 = afftile(:,:,idx_tm1);
    tstats_tm1 = pixstats(idx_tm1,:); % [idxt  idxtp1 med[idxt,idxtp1] min[idxt,idxtp1] max[idxt,idxtp1]];
end
if ~isnan(idx_tp1)
    tformtp1 = afftile(:,:,idx_tp1);
    tstats_tp1 = pixstats(idx_tp1,:); % [idxt  idxtp1 med[idxt,idxtp1] min[idxt,idxtp1] max[idxt,idxtp1]];
end

if isnan(idx_tp1) % to adj tiles below
    %% No tile below
    % assign a default value - defined in vectorField3D.m
    zlim_3 = zlimdefaults(3);%dims(3)-16;
    zlim_4 = zlimdefaults(4);%dims(3)-6;
    
    corr_ctrl_pnt_tmp(:,3) = zlim_3;
    % Position of the voxel points after affine transformation 
    contr_t_bot1 = [corr_ctrl_pnt_tmp, ones(Npts,1)] * tformt';
    % Compute the local displacement field by nearest neighbor
    % interpolation
    vect_bot1 = [Fxt(contr_t_bot1), Fyt(contr_t_bot1), Fzt(contr_t_bot1)];
    if noopt
        vect_bot1(:) = 0;
    end
    % Correct "physical" position of the control point
    contr_t_bot1_shifted = contr_t_bot1 + vect_bot1;
    
    corr_ctrl_pnt_tmp(:,3) = zlim_4;
    contr_t_bot2 = [corr_ctrl_pnt_tmp, ones(Npts,1)] * tformt';
    vect_bot2 = [Fxt(contr_t_bot2), Fyt(contr_t_bot2), Fzt(contr_t_bot2)];
    if noopt
        vect_bot2(:) = 0;
    end
    contr_t_bot2_shifted = contr_t_bot2 + vect_bot2;
    
    % book it
    %control(2*Npts+1:end,:,idxt) = [contr_t_bot1_shifted;contr_t_bot2_shifted];
    control_t_bot12 = [contr_t_bot1_shifted; contr_t_bot2_shifted];
else
    % Has a tile blow
    % check if there exists any descriptors between these two
    % tiles, if yes optimize to cover them, if not use a default
    % location
    % tstats: 1/ 2 / 3    / 4      / 5    / 6      / 7    / 8
    % tstats: t/tp1/med(t)/med(tp1)/min(t)/min(tp1)/max(t)/max(tp1)
    % get control lcocations
    if isnan(tstats(1)) % get a section from the bottom of tile as we have noisy samples at the top of a slice
        % will fill these after optimizing others
        nooptim = true;
        return
        %         % same as above for the bottom of the current tile t.
        %         zlim_3 = zlimdefaults(3);%dims(3)-11;
        %         zlim_4 = zlimdefaults(4);%dims(3)-1;
        %         tstats(8) = max(pixstats(:,8));
        %         tstats(6) = min(pixstats(:,6));
    else
        % bot1: then get top of overlap on layer tm1
        % (3)
        % Q: ???? (3) ?? meaning of the number
        zlim_3 = min(dims(3) - 2, tstats(5)); % min(tb.dims(idxt,3)-2, so that zlim_4 is bounded by zlim3<zlim4<=dims-1
        % bot2: then get bottom of overlap on layer tm1
        % (4)
        zlim_4 = min(dims(3) - 1, max(zlim_3 + 1, tstats(7))); % -2, as -1 results in error due to a bug in render ???
        if zlim_3 <= zlim_cntrl(2,idxt) % means that a tile is fully covered by the two adjacend tiles
            % @@ HEURISTIC
            zlim_cntrl(2,idxt) = max(1,zlim_3-1);
            zlim_cntrl(1,idxt) = max(0,zlim_3-6);
        end
    end
    
    corr_ctrl_pnt_tmp(:,3)= zlim_3;
    % Apply affine transformation
    contr_t_bot1 = [corr_ctrl_pnt_tmp ones(Npts,1)] * tformt';
    % Compute deformation field 
    vect_bot1=[Fxt(contr_t_bot1), Fyt(contr_t_bot1), Fzt(contr_t_bot1)];
    if noopt % Q: When is noopt = true? 
        vect_bot1(:) = 0;
    end
    contr_t_bot1_shifted = contr_t_bot1 + vect_bot1;
    
    corr_ctrl_pnt_tmp(:,3)= zlim_4;
    contr_t_bot2 = [corr_ctrl_pnt_tmp ones(Npts,1)] * tformt';
    vect_bot2 = [Fxt(contr_t_bot2), Fyt(contr_t_bot2), Fzt(contr_t_bot2)];
    if noopt
        vect_bot2(:) = 0;
    end
    contr_t_bot2_shifted = contr_t_bot2 + vect_bot2;
    % book it
    control_t_bot12 = [contr_t_bot1_shifted; contr_t_bot2_shifted];
    %% %%%%%%%%%%%%%%%%%%%%%
    % Search for the plane in the tile below that maximize the overlapping
    % region while making sure the corrected control point position is
    % above thecontr_t_bot2_shifted computed above
    for zlim_2 = (tstats(8)-1) : -1 : max(1,min(tstats(8)-1, tstats(6))) % scan from the maximum to the minumum, step -1
        % zlim_2 = max(zlim_1_bound + 1, tstats(8) - hbuf);
        corr_ctrl_pnt_tmp(:,3) = zlim_2;
        contr_tp1_top2 = [corr_ctrl_pnt_tmp, ones(Npts,1)]  * tformtp1';
        % shift it
        vectp1_top2 = [Fxtp1(contr_tp1_top2), Fytp1(contr_tp1_top2), Fztp1(contr_tp1_top2)];
        if noopt
            vectp1_top2(:) = 0;
        end
        contr_tp1_top2_shifted = contr_tp1_top2 + vectp1_top2;
        % make sure that there is overlap for all control points
        if all(contr_t_bot2_shifted(:,3) - contr_tp1_top2_shifted(:,3) > 0)
            %zlim_2
            break
        end
    end
    % top1: then get top of overlap on layer t at crop "h"
    % (1)
    % Q: No idea what is "h"...
    % htop is a preset parameter defined in main ( set to 5)
    zlim_1 = max(0, zlim_2 - htop);
    corr_ctrl_pnt_tmp(:,3) = zlim_1;
    contr_tp1_top1 = [corr_ctrl_pnt_tmp, ones(Npts,1)] * tformtp1';
    % shift it
    vectp1_top1 = [Fxtp1(contr_tp1_top1), Fytp1(contr_tp1_top1), Fztp1(contr_tp1_top1)];
    if noopt
        vectp1_top1(:) = 0;
    end
    contr_tp1_top1_shifted = contr_tp1_top1 + vectp1_top1;
    control_tp1_top12 = [contr_tp1_top1_shifted; contr_tp1_top2_shifted];
    zlim_cntrl(1:2,idx_tp1) = [zlim_1; zlim_2];
    zlim_cntrl(3:4,idxt) = [zlim_3; zlim_4];
end