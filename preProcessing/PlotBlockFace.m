%% Plot the bathymetry with blocked cell faces, ZhiLi20171021
function PlotBlockFace(bathCenter, block, surf, plotSurf, crange)
    dim = size(bathCenter);
    Dim = [size(block.Np,1) size(block.Np,2)];
    r = dim(1) / Dim(1);
    figure(1);
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Position',[100 100 800 400]);
    set(gcf,'Color',[1 1 1]);
    % extract the blocking face to plot
    bFaceNp = block.Np;
    bFaceOp = block.Op;
    bFaceNm = block.Nm;
    bFaceOm = block.Om;
    [val, ind] = min(abs(surf - plotSurf));
    bFaceNp = bFaceNp(:,:,ind);
    bFaceOp = bFaceOp(:,:,ind);
    bFaceNm = bFaceNm(:,:,ind);
    bFaceOm = bFaceOm(:,:,ind);
    % interpolate blocking faces on fine grids
    blockF = zeros(dim);
    for ii = 1:Dim(1)
        for jj = 1:Dim(2)
            if bFaceNp(ii,jj) == 1
                blockF(ii*r,(jj-1)*r+1:jj*r) = 10;                  
            end
            if bFaceOp(ii,jj) == 1
                blockF((ii-1)*r+1:ii*r,jj*r) = 10;
            end
            if bFaceNm(ii,jj) == 1
                blockF((ii-1)*r+1,(jj-1)*r+1:jj*r) = 10;                  
            end
            if bFaceOm(ii,jj) == 1
                blockF((ii-1)*r+1:ii*r,(jj-1)*r+1) = 10;
            end
        end
    end
    % overlap the blocking faces with the fine grid bathymetry
    aa = (blockF == 10);
    bathBlock = bathCenter;
    bathBlock(aa) = blockF(aa);
    % plot the fine grid bathymetry with blocking faces
    imagesc(bathBlock);
    caxis(crange);
    colormap(jet);
    colorbar;
    
    % plot the submerged bathymetry as comparison
    figure(2);
    set(gcf,'PaperPositionMode','auto');
    set(gcf,'PaperOrientation','landscape');
    set(gcf,'Position',[150 150 800 400]);
    set(gcf,'Color',[1 1 1]);
    bb = (bathCenter > plotSurf);
    bathCenter(bb) = 1;
    bathCenter(~bb) = 0;
    for ii = 1:r:dim(1)
        bathCenter(ii,:) = 0.2;
    end
    for jj = 1:r:dim(2)
        bathCenter(:,jj) = 0.2;
    end
    blockF(aa) = 0.7;
    bathCenter(aa) = blockF(aa);
    xVec = [1:dim(1)]/r;
    yVec = [1:dim(2)]/r;
    imagesc(yVec, xVec, bathCenter);
    colormap(jet);
end

