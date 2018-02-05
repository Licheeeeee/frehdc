%% Create the subgrid areas and volumes for the SubFREHD-C, ZhiLi20170709
% Added the automatic block-checking method by ZhiLi 20171027

%% Settings
% coarse and fine grid size, assume dx=dy
Dx = 10;
dx = 1;
r = Dx / dx;
% fine grid area
dA = dx^2;
% range of possible water elevations
surfmin = 0.7;
surfmax = 1.3;
dsurf = 0.005;
% mode for computing cell face areas, can be exact or min
faceMode = 'exact';
checkBlock = 0;
% whether or not use the Yh bottom drag model
useYh = 1;
YhReduction = 0;
dmax = 10;
% whether or not plot the face blocking for a specific face
plotBlock = 0;
plotSurf = 0.5;
crange = [-0.5 2.5];
% input fine bathymetry file name
fnameIn = 'bath_1x1_StairV3.mat';
% output file name
fnameA = 'subArea_1x1to10x10StairV3.mat';
fnameB = 'subBath_1x1to10x10StairV3.mat';

%% Load files
load(fnameIn);
% dimension for the fine and coarse bathymetry
dim = size(bathCenter);
Dim = dim / r;
% vector of all possible water elevations
surf = surfmin:dsurf:surfmax;
N = length(surf);
% curvature
CvX = zeros(dim);
CvY = zeros(dim);

%% Create all required subgrid fields
subA.surf = surf';
subA.dx = dx;
subA.Dx = Dx;
subA.V = zeros([Dim N]);
subA.Z = zeros([Dim N]);
subA.N = zeros([Dim N]);
subA.O = zeros([Dim N]);
subA.Vxp = zeros([Dim N]);
subA.Vyp = zeros([Dim N]);
subA.Vxm = zeros([Dim N]);
subA.Vym = zeros([Dim N]);
subA.Zxp = zeros([Dim N]);
subA.Zyp = zeros([Dim N]);
subA.Zxm = zeros([Dim N]);
subA.Zym = zeros([Dim N]);
subA.Np = zeros([Dim N]);
subA.Op = zeros([Dim N]);
subA.Nm = zeros([Dim N]);
subA.Om = zeros([Dim N]);
subA.CvX = zeros([Dim N]);
subA.CvY = zeros([Dim N]);
subA.Yh = zeros([Dim N]);
subB.bottom = zeros(Dim);
subB.bottomXP = zeros(Dim);
subB.bottomYP = zeros(Dim);
subB.wdOp = zeros(Dim);
subB.wdOm = zeros(Dim);
subB.wdNp = zeros(Dim);
subB.wdNm = zeros(Dim);

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
    [block] = CheckFaceBlocking(bathCenter, surf, Dx/dx);
end

%% Compute the subgrid areas and volumes for each possible elevation
for kk = 1:N
    h = surf(kk);
    fprintf('Computing subgrid variables for elevation = %f...\n',h);
    depth = max(h - bathCenter, 0);
    % compute the curvature only for the wet regions
    if r > 1
        CvXfine = zeros(dim);
        CvYfine = zeros(dim);
        depthCntr = depth(2:end-1,2:end-1);
        depthXp = depth(3:end,2:end-1);
        depthXm = depth(1:end-2,2:end-1);
        depthYp = depth(2:end-1,3:end);
        depthYm = depth(2:end-1,1:end-2);
        isdry = (depthCntr <= 0);   depthCntr(isdry) = NaN;
        isdry = (depthXp <= 0);   depthXp(isdry) = NaN;
        isdry = (depthXm <= 0);   depthXm(isdry) = NaN;
        isdry = (depthYp <= 0);   depthYp(isdry) = NaN;
        isdry = (depthYm <= 0);   depthYm(isdry) = NaN;

        CvXfine(2:end-1,2:end-1) = ...
            abs((depthXp - 2*depthCntr + depthXm) / (2*dx));
        CvYfine(2:end-1,2:end-1) = ...
            abs((depthYp - 2*depthCntr + depthYm) / (2*dx));
        aa = isnan(CvXfine);   CvXfine(aa) = 0;
        aa = isnan(CvYfine);   CvYfine(aa) = 0;
    end
    
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
            % compute coarse grid curvatures
            if r > 1
                r2 = round(r/2);
                if ii ~= Dim(1)
                    localCvX = CvXfine(i1+r2:i2+r2,j1:j2);
                    subA.CvX(ii,jj,kk) = sum(localCvX(:));
                end
                if jj ~= Dim(2)
                    localCvY = CvYfine(i1:i2,j1+r2:j2+r2);
                    subA.CvY(ii,jj,kk) = sum(localCvY(:));
                end
            end
            
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
            % compute cell center areas
            subA.N(ii,jj,kk) = dx .* min(sum(localDepth, 2));
            subA.O(ii,jj,kk) = dx .* min(sum(localDepth, 1));
%             subA.N(ii,jj,kk) = dx .* sum(localDepth(r/2,:));
%             subA.O(ii,jj,kk) = dx .* sum(localDepth(:,r/2));
            if r > 1
                % compute cell face volumes
                subA.Vxp(ii,jj,kk) = dA .* sum(sum(localDepth(round(r/2)+1:r,:)));
                subA.Vyp(ii,jj,kk) = dA .* sum(sum(localDepth(:,round(r/2)+1:r)));
                subA.Vxm(ii,jj,kk) = dA .* sum(sum(localDepth(1:round(r/2),:)));
                subA.Vym(ii,jj,kk) = dA .* sum(sum(localDepth(:,1:round(r/2))));
                % compute cell face free surface area
                subA.Zxp(ii,jj,kk) = dA .* sum(sum(localDepth(round(r/2)+1:r,:) > 0));
                subA.Zyp(ii,jj,kk) = dA .* sum(sum(localDepth(:,round(r/2)+1:r) > 0));
                subA.Zxm(ii,jj,kk) = dA .* sum(sum(localDepth(1:round(r/2),:) > 0));
                subA.Zym(ii,jj,kk) = dA .* sum(sum(localDepth(:,1:round(r/2)) > 0));
%                 % check for blocked cell faces
                if checkBlock == 1 && r > 1
                    % remove face blockings
                    if block.Np(ii,jj,kk) == 1
                        if subA.Np(ii,jj,kk) == 0
%                                 block.Np(ii,jj,kk) = 0;
                        else
                            subA.Np(ii,jj,kk) = 0;
                        end
                    end
                    if block.Nm(ii,jj,kk) == 1
                        if subA.Nm(ii,jj,kk) == 0
%                                 block.Nm(ii,jj,kk) = 0;
                        else
                            subA.Nm(ii,jj,kk) = 0;
                        end
                    end
                    if block.Op(ii,jj,kk) == 1
                        if subA.Op(ii,jj,kk) == 0
%                                 block.Op(ii,jj,kk) = 0;
                        else
                            subA.Op(ii,jj,kk) = 0;
                        end
                    end
                    if block.Om(ii,jj,kk) == 1
                        if subA.Om(ii,jj,kk) == 0
%                                 block.Om(ii,jj,kk) = 0;
                        else
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
                % compute cell face free surface area
                subA.Zxp(ii,jj,kk) = dA .* (localDepth > 0) ./ 2;
                subA.Zyp(ii,jj,kk) = dA .* (localDepth > 0) ./ 2;
                subA.Zxm(ii,jj,kk) = dA .* (localDepth > 0) ./ 2;
                subA.Zym(ii,jj,kk) = dA .* (localDepth > 0) ./ 2;
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
%                 aa = (localDepth > 0);
%                 M = sum(aa(:));
%                 local3 = localDepth(aa).^(1/3);
%                 if M > 0
%                     subA.Yh(ii,jj,kk) = nansum(1 ./ local3) ./ M;
%                 end
                
                
            end
        end
    end
    
    % reduce bottom drag due to insufficient grid resolutions in channels
    if YhReduction == 1
        [subA.Yh] = SubgridDragReduction(subA.Yh, subA.V, dmax, kk);
    end
    
end


%% Save output

save(fnameA, 'subA' ,'-v7.3');
save(fnameB, 'subB', '-v7.3');

%% Plot the blocking face
if plotBlock == 1
    PlotBlockFace(bathCenter, block, surf, plotSurf, crange);
end









