function [vecfield] = vectorField3D(params,scopeloc,regpts,scopeparams,curvemodel,theselayers)
%VECTORFIELD3D Summary of this function goes here
%
% [OUTPUTARGS] = VECTORFIELD3D(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/26 19:53:38 $	$Revision: 0.1 $
% Copyright: HHMI 2017
%%
disp('Calculating vector fields');
if nargin<5
    affcase = 2;
end
noopt = 0;
htop = params.htop;
Nlayer = params.Nlayer;
Npts = (params.Ndivs+1).^2;
dims = params.imagesize;
order = params.order;
zlimdefaults = [2 20 dims(3)-21 dims(3)-3];
params.zlimdefaults = zlimdefaults;
beadparams = params.beadparams;
params.applyFC=1;
%%
[tileneighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% tileneighbors = tileneighbors(:,checkthese);
%% % initialize points
numTiles = length(scopeloc.filepath);
% Determine the starting and ending pixel subscript for generating the
% control points
[st,ed] = util.getcontolpixlocations(scopeloc,params,scopeparams);
% Generate control points (in a regular grid) 
[corrctrlpnttmp_, xlim_cntrl, ylim_cntrl] = util.getCtrlPts(params.imagesize(1),params.imagesize(2),params,st,ed);
% zlimdefaults = [5 25 dims(3)-26 dims(3)-6];
subcorrctrlpnttmp = corrctrlpnttmp_ + 1;
idxctrl = sub2ind(dims([2 1]), subcorrctrlpnttmp(:,2), subcorrctrlpnttmp(:,1));
%% Use the affine matrix after field curvature correction 
afftransform = reshape([scopeparams(:).affineglFC],3,3,[]);
% flip x/y dimensions to match imaging
mflip = -eye(3);%mflip(3,3)=1;
for ii=1:size(afftransform,3)
    afftransform(:,:,ii) = afftransform(:,:,ii)*mflip;
end

%% for efficiency reasons precompute this
xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];

%%
control = zeros(Npts * Nlayer, 3, numTiles);
top_cntrl_1 = zeros(Npts, 3, numTiles);
top_cntrl_2 = zeros(Npts, 3, numTiles);
bot_cntrl_1 = zeros(Npts, 3, numTiles);
bot_cntrl_2 = zeros(Npts, 3, numTiles);
zlim_cntrl = zeros(Nlayer, numTiles);
afftile = zeros(3, 4, numTiles);
pixstats = nan(numTiles,8);

filep = strsplit(scopeloc.filepath{1},'/');
vecfield.root = fullfile('/',filep{1:end-4}); % The root is set according to the scopleloc 

%% Initialization 
for idxt = 1:numTiles

    % copy path
    filepath = fileparts(scopeloc.filepath{idxt});
    vecfield.path{idxt} = filepath((length(vecfield.root)+1) : end);
    % form an affine matrix - the last column is for translation 
    tform = [afftransform(:,:,idxt) scopeloc.loc(idxt,:)'*1e6];
    afftile(:,:,idxt) = tform;
    
    % field curvature
    [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,idxt),order,xy,dims,subcorrctrlpnttmp);
    corr_ctrl_pnt_tmp = locs-1;
    
    % top layer at z=15: as we have some noisy data at the top
    zlim_1 = zlimdefaults(1);
    corr_ctrl_pnt_tmp(:,3)= zlim_1;
    contr_t_top1 = [corr_ctrl_pnt_tmp ones(Npts,1)] * tform';
    
    % second layer at z=25
    zlim_2 = zlimdefaults(2);
    corr_ctrl_pnt_tmp(:,3)= zlim_2;
    contr_t_top2 = [corr_ctrl_pnt_tmp ones(Npts,1)] * tform';
    
    % third layer at z=251-26
    zlim_3 = zlimdefaults(3);
    corr_ctrl_pnt_tmp(:,3)= zlim_3;
    contr_t_bot1 = [corr_ctrl_pnt_tmp ones(Npts,1)] * tform';
    
    % forth layer at z=251-16
    zlim_4 = zlimdefaults(4);
    corr_ctrl_pnt_tmp(:,3)= zlim_4;
    contr_t_bot2 = [corr_ctrl_pnt_tmp ones(Npts,1)] * tform';
    
    top_cntrl_1(:,:,idxt) = contr_t_top1;
    top_cntrl_2(:,:,idxt) = contr_t_top2;
    bot_cntrl_1(:,:,idxt) = contr_t_bot1;
    bot_cntrl_2(:,:,idxt) = contr_t_bot2;
    
    control(1:2*Npts,:,idxt) = [contr_t_top1; contr_t_top2];
    control(2*Npts+1:end,:,idxt) = [contr_t_bot1; contr_t_bot2];
    zlim_cntrl(:,idxt) = [zlim_1 zlim_2 zlim_3 zlim_4];
    
    % pix stats
    layer = regpts{idxt}.X;
    layerp1 = regpts{idxt}.Y;
    if isempty(layer)
        continue
    end
    meddesc = round(median([layer(:,3) layerp1(:,3)],1));% median z-coordinate in tile 1 and tile 2
    mindesc = round(min([layer(:,3) layerp1(:,3)],[],1));% min z-coordinate in tile 1 and tile 2
    maxdesc = round(max([layer(:,3) layerp1(:,3)],[],1));% max z-coordinate in tile 1 and tile 2
    pixstats(idxt,:) = [idxt  regpts{idxt}.neigs(4) meddesc mindesc maxdesc];
    % pixstats(62,2:end) = pixstats(1162,2:end);
end
%%
if 1
    Rmin = min(scopeloc.loc)*1e6; 
    Rmin = Rmin * 0.995; % Origin of the rendering space - expand by 0.005
    Rmax = round((max(scopeloc.loc) * 1e3 + scopeparams(1).imsize_um) * 1e3);
    Rmax = Rmax*1.005; % Maximum coordinate values of the rendering space - expand by 0.005
    targetidx = [];
    latticeZRange = unique(scopeloc.gridix(:,3));
    if nargin<6 || isempty(theselayers)
        theselayers = latticeZRange(1:end-1)';
    end
%     stopindex = 149; % Q: meaning?
    nooptim = NaN(1,numTiles);
    for iter_layer = theselayers
        %%
        disp(['    Layer ' num2str(iter_layer) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
        ix = (scopeloc.gridix(:,3)'==iter_layer);
        idx_in_layer = find(ix);
        if (sum(ix)<1)
            disp(['No tiles found in layer #' num2str(iter_layer) '!!']);
            continue;
        end
        %% get interpolants based on paired descriptors
        [Fxt, Fyt, Fzt, Fxtp1, Fytp1, Fztp1,XYZ_tori,XYZ_tp1ori,outliers] =...
            util.getInterpolants(ix, regpts, afftile, params, curvemodel);
        %% Visualization 
        for dum2 = 1
            if params.debug | 0
                %%
                [bandwidth,density,X,Y] = kde2d(XYZ_tori(:,1:2),256,Rmin(1:2),Rmax(1:2));
                figure(35), clf, cla,
                view([0,90])
                set(gca,'Ydir','reverse')
                hold on
                colormap hot,
                set(gca, 'color', 'w');
                %         plot(XYZ_tori(:,1),XYZ_tori(:,2),'w+','MarkerSize',5)
                
                plot3(XYZ_tori(:,1),XYZ_tori(:,2),XYZ_tori(:,3),'k.') % layer t
                plot3(XYZ_tp1ori(:,1),XYZ_tp1ori(:,2),XYZ_tp1ori(:,3),'r.') % layer tp1
                myplot3(XYZ_tori(outliers,:),'md') % layer t
                myplot3(XYZ_tori(outliers,:),'gd') % layer t
                x = scopeloc.loc(ix,1)*1e6;
                y = scopeloc.loc(ix,2)*1e6;
                w = scopeparams(1).imsize_um(1)*ones(sum(ix),1)*1e3;
                h = scopeparams(1).imsize_um(2)*ones(sum(ix),1)*1e3;
                for ii=1:sum(ix)
                    rectangle('Position', [x(ii)-w(ii) y(ii)-h(ii) w(ii) h(ii)],'EdgeColor','r')
                    text(x(ii)-w(ii)/2,y(ii)-h(ii)/2,2,num2str(ii),...
                        'Color','m','HorizontalAlignment','center','FontSize',8)
                end
                
                for ii=find(idx_in_layer==5163)
                    rectangle('Position', [x(ii)-w(ii) y(ii)-h(ii) w(ii) h(ii)],'EdgeColor','w','LineWidth',2,'FaceColor',[0 .5 .5])
                end
                %surfc(X,Y,density/max(density(:)),'LineStyle','none'),
                set(gca,'Ydir','reverse')
                %         legend('P_t','P_{t+1}','E_t','E_{t+1}')
                set(gca, ...
                    'Box'         , 'on'     , ...
                    'TickDir'     , 'out'     , ...
                    'TickLength'  , [.02 .02]/2 , ...
                    'XTickLabel'  , [1:10] , ...
                    'YTickLabel'  , [1:10] , ...
                    'XMinorTick'  , 'on'      , ...
                    'YMinorTick'  , 'on'      , ...
                    'XGrid'       , 'on'      , ...
                    'YGrid'       , 'on'      , ...
                    'XColor'      , [.3 .3 .3], ...
                    'YColor'      , [.3 .3 .3], ...
                    'LineWidth'   , 1         );
                ax = gca;
                ax.YLabel.String = 'mm';
                ax.YLabel.FontSize = 16;
                ax.XLabel.String = 'mm';
                ax.XLabel.FontSize = 16;
                xlim([Rmin(1) Rmax(1)]-params.imsize_um(1)*1e3)
                ylim([Rmin(2) Rmax(2)]-params.imsize_um(2)*1e3)
                % title([num2str(t),' - '])
                drawnow
                %%
                vizfolder = 'qualityfold'
                if ~exist(vizfolder,'dir')
                    mkdir(vizfolder)
                end
                export_fig(fullfile(vizfolder,sprintf('Slice-%05d.png',iter_layer)))
            end
        end
        %% Print information 
        if isempty(Fxt) | size(Fxt.Points,1)<10,
            disp(['    MISSING SLICE @ Layer ' num2str(iter_layer) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
            continue,
        else
            disp(['    Layer ' num2str(iter_layer) ' of ' num2str(max(scopeloc.gridix(:,3))) ' totDesc: ' num2str(size(Fxt.Points,1))]);
        end
        %%
%         if any(idx_in_layer==stopindex)
%             stop=1;
%         end
        %% Compute the position of the corresponding control points in two tiles
        for idxt = idx_in_layer % layer t
%             if any(idxt==stopindex)
%                 stop=1;
%             end
            if ~isempty(targetidx)
                if ~any(idxt==targetidx)
                    continue
                end
            end
            idxtm1 = tileneighbors(idxt,6);
            idxtp1 = tileneighbors(idxt,7);
            
            [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,idxt),order,xy,dims,subcorrctrlpnttmp);
            corr_ctrl_pnt_tmp = locs-1;
            [nooptim(idxt),control_t_bot12,control_tp1_top12,zlim_cntrl] = optimpertile...
                (idxt,params,tileneighbors,afftile,pixstats,zlim_cntrl,corr_ctrl_pnt_tmp,Fxt,Fyt,Fzt,Fxtp1,Fytp1,Fztp1);
            
            if nooptim(idxt) % nomatch with below - in which case can this happen ? 

            elseif isempty(control_tp1_top12) % no below adjacent tile
                control(2*Npts+1:end,:,idxt) = control_t_bot12;
            else
                control(2*Npts+1:end,:,idxt) = control_t_bot12;
                control(1:2*Npts,:,idxtp1) = control_tp1_top12;
            end
            
        end
        %% Deal with nooptim 
        % KNN search to get / replace the missing parameters
        anchorinds = idx_in_layer(~nooptim(idx_in_layer));
        anchors = scopeloc.gridix(anchorinds,:);        
        queryinds = idx_in_layer(nooptim(idx_in_layer)>0);
        queries = scopeloc.gridix(queryinds,:);
        IDX_nn = knnsearch(anchors,queries,'K',1);
        [IDX,dist] = rangesearch(anchors,queries,sqrt(2));
        for ine = 1 : length(queryinds)
            if isempty(IDX{ine})
                IDX{ine} = IDX_nn(ine);
            end
        end        
        %
        for idxt = idx_in_layer % layer t
            idxtm1 = tileneighbors(idxt,6);
            idxtp1 = tileneighbors(idxt,7);
            if nooptim(idxt)
                indquer = find(idxt==queryinds);
                zlimdefaults = round(median(zlim_cntrl(:,anchorinds(IDX{indquer}) ),2))';
                
                [nooptim(idxt),control_t_bot12,control_tp1_top12,zlim_cntrl] = nomatchoptim...
                    (idxt,params,tileneighbors,afftile,pixstats,zlim_cntrl,corr_ctrl_pnt_tmp,...
                    Fxt,Fyt,Fzt,Fxtp1,Fytp1,Fztp1,zlimdefaults);

                control(2*Npts+1:end,:,idxt) = control_t_bot12;
                control(1:2*Npts,:,idxtp1) = control_tp1_top12;
            end
        end
    end
    bbox = zeros(8,3,numTiles);
    origin = zeros(numTiles,3);
    sz = zeros(numTiles,3);
    for i = 1:numTiles
        [bbox(:,:,i), origin(i,:), sz(i,:)] = util.bboxFromCorners(control(:,:,i));
    end
%% update affines
    disp('Update affines based on control points')
    tform = zeros(5,5,numTiles);
    numX = size(control(:,:,1),1);
    for ii = 1:numTiles
        X = control(:,:,ii)'/1000;
        x = xlim_cntrl;
        y = ylim_cntrl;
        z = zlim_cntrl(:,ii)';
        [xx,yy,zz] = ndgrid(x,y,z);
        x = [xx(:),yy(:),zz(:)]';
        X = [X;ones(1,numX)];
        x = [x;ones(1,numX)];
        %A = (X*x')/(x*x');
        A = X/x;
        Aest = eye(5);
        Aest(1:3,1:3) = A(1:3,1:3)*1000;
        Aest(1:3,5) = A(1:3,4)*1000;
        Aest(5,1:3) = A(4,1:3)*1000;
        tform(:,:,ii) = Aest;
    end
%%  Output
    vecfield.control = control;
    vecfield.xlim_cntrl = xlim_cntrl;
    vecfield.ylim_cntrl = ylim_cntrl;
    vecfield.zlim_cntrl = zlim_cntrl;
    vecfield.afftile = afftile;
    vecfield.tform = tform;
    vecfield.corrctrlpnttmp = corr_ctrl_pnt_tmp;
    vecfield.ctrlpnttmp = subcorrctrlpnttmp-1;
    vecfield.bbox = bbox;
    vecfield.origin = origin;
    vecfield.sz = sz;
    vecfield.time = datestr(now);
    vecfield.theselayers = theselayers;
    
    
end














