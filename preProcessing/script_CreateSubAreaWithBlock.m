%% Create the subgrid areas and volumes for the SubFREHD-C, ZhiLi20170709
%  ###Version1: Added the automatic block-checking method by ZhiLi 20171027
%  ###Version2: The small pixels in block-checking function are removed by
%  ZhiLi20180213. 
%  ###Version3: Added an effective sidewall drag model by ZhiLi 20180327
%  ###Version4: Added an effective bottom drag model by ZhiLi 20180409
%  ###Version7: Compute slope rather than curvature for bottom friction
%  correction. Version5 and 6 are abandoned. ZhiLi 20180422
%  ###Version8: Added the Manning's coefficient by Casas2010. ZhiLi20180506
%  ###Version9: Enabled subgrid variables only for part of the domain.
%  ZhiLi20180621
%% Settings
% coarse and fine grid size, assume dx=dy
Dx = 20;
dx = 1;
r = Dx / dx;
% fine grid area
dA = dx^2;
% range of possible water elevations
surfmin = -0.2;
surfmax = 0.8;
dsurf = 0.01;
% mode for computing cell face areas, can be exact or min
faceMode = 'exact';
checkBlock = 1;
% whether or not use the Yh bottom drag model
useYh = 0;
% whether or not use effective sidewall drag Cd correction
useEffCd = 1;
% whether or not use Casas model for Manning's n
useCn = 0;
n0 = 0.03;
% input fine bathymetry file name
% fnameIn = 'nuecesUpChannelBay_1x1.mat';
fnameIn = '../../nueces2007_frehdc1x1_finalV4_20180624.mat';
load(fnameIn);
% dimension for the fine and coarse bathymetry
dim = size(bathCenter);
Dim = dim / r;
% whether or not use subgrid method only for part of the domain
usePartialSub = 1;
pdomain = [1 245 1 330];
fnameCoarse = 'NDHMV4_Bath_Edges_and_Channel_Final_dH0p2_filterSize_40x40_gridSize_20x20.mat';
% whether or not plot the face blocking for a specific face
plotBlock = 0;
plotSurf = 0.5;
crange = [-0.5 2.5];
% output file name
fnameA = 'subArea_NDHM20x20V4FEC.mat';
fnameB = 'subBath_NDHM20x20V4FEC.mat';

%% Load files
if usePartialSub == 1
    bathCoar = load(fnameCoarse);
end
% vector of all possible water elevations
surf = surfmin:dsurf:surfmax;
N = length(surf);

%% Create all required subgrid fields
subA.surf = surf';
subA.dx = dx;
subA.Dx = Dx;
subA.V = zeros([Dim N]);
subA.Z = zeros([Dim N]);
subA.Vxp = zeros([Dim N]);
subA.Vyp = zeros([Dim N]);
subA.Vxm = zeros([Dim N]);
subA.Vym = zeros([Dim N]);
subA.Np = zeros([Dim N]);
subA.Op = zeros([Dim N]);
subA.Nm = zeros([Dim N]);
subA.Om = zeros([Dim N]);
subA.CvX = zeros([Dim N]);
subA.CvY = zeros([Dim N]);
subA.effCdX = zeros([Dim N]);
subA.effCdY = zeros([Dim N]);
subB.bottom = zeros(Dim);
subB.bottomXP = zeros(Dim);
subB.bottomYP = zeros(Dim);
subB.wdOp = zeros(Dim);
subB.wdOm = zeros(Dim);
subB.wdNp = zeros(Dim);
subB.wdNm = zeros(Dim);
if useYh == 1
    subA.Yh = zeros([Dim N]);
end
if useCn == 1
    subA.Cn = zeros([Dim N]);
end

%% Compute the subgrid bathymetry
for ii = 1:Dim(1)
    for jj = 1:Dim(2)
        i1 = (ii-1)*r+1;
        i2 = ii*r;
        j1 = (jj-1)*r+1;
        j2 = jj*r;
        bath = bathCenter(i1:i2,j1:j2);
        % bottom elevation of the coarse grid
        subB.bottom(ii,jj) = min(bath(:));
   
        if i2 == dim(1)
            localBathx = bathCenter(i2, j1:j2);
        else
            localBathx = bathCenter(i2:i2+1, j1:j2);
            localBathP = min(localBathx(1,:));
            localBathM = min(localBathx(2,:));
        end
        if j2 == dim(2)
            localBathy = bathCenter(i1:i2,j2);
        else
            localBathy = bathCenter(i1:i2,j2:j2+1);
            localBathP = min(localBathy(:,1));
            localBathM = min(localBathy(:,2));
        end
        
        subB.bottomXP(ii,jj) = max(min(localBathx,[],2));
        subB.bottomYP(ii,jj) = max(min(localBathy,[],1));
        % compute the wetting drying elevation
        if r > 1
            halfpiece = bath(round(r/2)+1:r,:);
        else
            halfpiece = bath;
        end
        [val, ind] = max(sum(halfpiece, 2));
        subB.wdNp(ii,jj) = min(halfpiece(ind,:));
        
        if r > 1
            halfpiece = bath(1:round(r/2),:);
        else
            halfpiece = bath;
        end
        [val, ind] = max(sum(halfpiece, 2));
        subB.wdNm(ii,jj) = min(halfpiece(ind,:));
        
        if r > 1
            halfpiece = bath(:,round(r/2)+1:r);
        else
            halfpiece = bath;
        end
        [val, ind] = max(sum(halfpiece, 1));
        subB.wdOp(ii,jj) = min(halfpiece(:,ind));
        
        if r > 1
            halfpiece = bath(:,1:round(r/2));
        else
            halfpiece = bath;
        end
        [val, ind] = max(sum(halfpiece, 1));
        subB.wdOm(ii,jj) = min(halfpiece(:,ind));
    end
end

%% Check for face blockings
if checkBlock == 1 && r > 1
    % compute effective side drag term, ZhiLi 20180327
    if useEffCd == 1
        [cd] = ComputeDragCorrection(bathCenter, surf, r);
        [cv] = ComputeBottomSlope(bathCenter, surf, r);
        subA.effCdX = cd.effCdX;
        subA.effCdY = cd.effCdY;
        subA.CvX = cv.CvX;
        subA.CvY = cv.CvY;
        clear cd cv
    end
    [block] = CheckFaceBlocking(bathCenter, surf, Dx/dx);
end

%% Compute the subgrid areas and volumes for each possible elevation
for kk = 1:N
    h = surf(kk);
    fprintf('Computing subgrid variables for elevation = %f...\n',h);
    depth = max(h - bathCenter, 0);
    % compute the cell face area
    if r > 1
        [subA] = ComputeSubgridFaceArea(subA, depth, dim, dx, r, kk);
    end
    
    % compute the variables on the cell faces
    for ii = 1:Dim(1)
        for jj = 1:Dim(2)
            i1 = (ii-1)*r+1;
            i2 = ii*r;
            j1 = (jj-1)*r+1;
            j2 = jj*r;
            localDepth = depth(i1:i2,j1:j2);
            % adjust at boundaries
            if i2 == dim(1)
                localDepthx = depth(i1+round(r/2):i2, j1:j2);
            else
                localDepthx = depth(i1+round(r/2):i2+round(r/2), j1:j2);
            end
            if j2 == dim(2)
                localDepthy = depth(i1:i2, j1+round(r/2):j2);
            else
                localDepthy = depth(i1:i2, j1+round(r/2):j2+round(r/2));
            end
            % compute cell center volumes
            subA.V(ii,jj,kk) = dA .* sum(localDepth(:));
            % compute cell center free surface area
            subA.Z(ii,jj,kk) = dA .* sum(sum(localDepth > 0));
            if r > 1
                % compute cell face volumes
                subA.Vxp(ii,jj,kk) = dA .* sum(sum(localDepth(round(r/2)+1:r,:)));
                subA.Vyp(ii,jj,kk) = dA .* sum(sum(localDepth(:,round(r/2)+1:r)));
                subA.Vxm(ii,jj,kk) = dA .* sum(sum(localDepth(1:round(r/2),:)));
                subA.Vym(ii,jj,kk) = dA .* sum(sum(localDepth(:,1:round(r/2))));
                % check for blocked cell faces
                if checkBlock == 1 && r > 1
                    % remove face blockings
                    if block.Np(ii,jj,kk) == 1
                        if subA.Np(ii,jj,kk) ~= 0
                            subA.Np(ii,jj,kk) = 0;
                        end
                    end
                    if block.Nm(ii,jj,kk) == 1
                        if subA.Nm(ii,jj,kk) ~= 0
                            subA.Nm(ii,jj,kk) = 0;
                        end
                    end
                    if block.Op(ii,jj,kk) == 1
                        if subA.Op(ii,jj,kk) ~= 0
                            subA.Op(ii,jj,kk) = 0;
                        end
                    end
                    if block.Om(ii,jj,kk) == 1
                        if subA.Om(ii,jj,kk) ~= 0
                            subA.Om(ii,jj,kk) = 0;
                        end
                    end
                    subA.block = block;
                end
            else
                subA.Vxp(ii,jj,kk) = dA .* localDepth ./ 2;
                subA.Vyp(ii,jj,kk) = dA .* localDepth ./ 2;
                subA.Vxm(ii,jj,kk) = dA .* localDepth ./ 2;
                subA.Vym(ii,jj,kk) = dA .* localDepth ./ 2;
                % compute cell face areas
                subA.Np(ii,jj,kk) = dx .* localDepth;
                subA.Op(ii,jj,kk) = dx .* localDepth;
                subA.Nm(ii,jj,kk) = dx .* localDepth;
                subA.Om(ii,jj,kk) = dx .* localDepth;
            end
            % compute Yh for the drag model, where Cd = Yh * gn^2
            if useYh == 1
                Y = subA.V(ii,jj,kk) ./ (Dx^2);
                aa = (localDepth > 0);
                if sum(aa(:)) > 0
                    sumh = (nansum(localDepth(aa) .^ (5/3))) .^ (-2);
                    subA.Yh(ii,jj,kk) = (Dx^2/dx^2) .* (Y^3) .* sumh;
                else
                    subA.Yh(ii,jj,kk) = 0;
                end
            end
            % compute Cn for Casas2010 model
            if useCn == 1
                D = max(0, mean(mean(bathCenter(i1:i2,j1:j2) - subB.bottom(ii,jj))));
                H = max(0, surf(kk) - subB.bottom(ii,jj));
                if D > 0
                    xi = H / D;
                else
                    xi = 100;
                end
                % check if xi is within the range
                if xi >= 0.2 && xi <= 7
                    fxi = 1 + (1/xi)*log(cosh(1-xi)/cosh(1));
                    subA.Cn(ii,jj,kk) = H^(1/6) / sqrt(9.81*4.5*fxi);
                else
                    subA.Cn(ii,jj,kk) = n0;
                end
            end
        end
    end
    
end


%% Remove subgrid variable out of the partial domain, ZhiLi20180621
if usePartialSub == 1
    for kk = 1:N
        h = surf(kk);
        fprintf('Disabling subgrid variables for h = %f...\n',h);
        for ii = 1:Dim(1)
            for jj = 1:Dim(2)
                if ii < pdomain(1) || ii > pdomain(2) ||...
                        jj < pdomain(3) || jj > pdomain(4)
                    b = bathCoar.bathCenter(ii,jj);
                    depth = max(h - b, 0);
                    subA.V(ii,jj,kk) = depth * Dx * Dx;
                    subA.Z(ii,jj,kk) = Dx * Dx;
                    subA.Vxp(ii,jj,kk) = 0.5 * subA.V(ii,jj,kk);
                    subA.Vyp(ii,jj,kk) = 0.5 * subA.V(ii,jj,kk);
                    subA.Vxm(ii,jj,kk) = 0.5 * subA.V(ii,jj,kk);
                    subA.Vym(ii,jj,kk) = 0.5 * subA.V(ii,jj,kk);
                    subA.Np(ii,jj,kk) = depth * Dx;
                    subA.Nm(ii,jj,kk) = depth * Dx;
                    subA.Op(ii,jj,kk) = depth * Dx;
                    subA.Om(ii,jj,kk) = depth * Dx;
                    subA.effCdX(ii,jj,kk) = 0;
                    subA.effCdY(ii,jj,kk) = 0;
                    subA.CvX(ii,jj,kk) = 0;
                    subA.CvY(ii,jj,kk) = 0;
                    if kk == 1
                        subB.bottom(ii,jj) = b;
                        if ii < Dim(1)
                            subB.bottomXP(ii,jj) = max(b, bathCoar.bathCenter(ii+1,jj));
                        else
                            subB.bottomXP(ii,jj) = b;
                        end
                        subB.wdNp(ii,jj) = subB.bottomXP(ii,jj);
                        if ii > 1
                            subB.wdNm(ii,jj) = subB.bottomXP(ii-1,jj);
                        else
                            subB.wdNm(ii,jj) = subB.bottomXP(ii,jj);
                        end
                        if jj < Dim(2)
                            subB.bottomYP(ii,jj) = max(b, bathCoar.bathCenter(ii,jj+1));
                        else
                            subB.bottomYP(ii,jj) = b;
                        end
                        subB.wdOp(ii,jj) = subB.bottomYP(ii,jj);
                        if jj > 1
                            subB.wdOm(ii,jj) = subB.bottomYP(ii,jj-1);
                        else
                            subB.wdOm(ii,jj) = subB.bottomYP(ii,jj);
                        end
                    end
                end
            end
        end
    end
end

%% Save output

save(fnameA, 'subA' ,'-v7.3');
save(fnameB, 'subB', '-v7.3');

%% Plot the blocking face
if plotBlock == 1
    PlotBlockFace(bathCenter, block, surf, plotSurf, crange);
end









