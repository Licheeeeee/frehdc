%% Check blocked cell faces, ZhiLi20171027
function [block] = CheckFaceBlocking(bathCenter, surf, r)
    % ---------- User Settings ----------
    % bathCenter = fine grid data
    % surf = all possible free surface elevations
    % dx = the grid size hierachy
    %      the first element is the input fine grid size
    %      the last element is the output coarse grid size
    % whether or not check inter-cell and cross-cell connectivity
    intercell = 1;
    crosscell = 1;
    % the input data size
    dim = size(bathCenter);
    % the output data size
    Dim = dim / r;    
    N = length(surf);
    % initialize output structure
    block.Np = zeros([Dim N]);
    block.Nm = zeros([Dim N]);
    block.Op = zeros([Dim N]);
    block.Om = zeros([Dim N]);
    % ------------ Computation ------------
    % loop over all possible surface elevations
    for kk = 1:N
        fprintf('Surface elevation = %f ...\n',surf(kk));        
        aa = (surf(kk) > bathCenter);
        bathBinary = zeros(dim);
        % wet area is 1, dry area is 0
        bathBinary(aa) = 1;
        bathBinary(~aa) = 0;
        % check for cell blockings
        for ii = 1:Dim(1)
            for jj = 1:Dim(2)
                % extract the binary data for each cell
                gridcell = bathBinary((ii-1)*r+1:ii*r,(jj-1)*r+1:jj*r);
                % ------ step 1: delete small wet areas ------
                % reduce the number of wet areas, we need only 1 wet area
                % in a coarse grid cell
                [bwater, Nw] = bwlabel(gridcell, 4);
                if Nw > 0
                    % find the largest wet cluster
                    waterSize = zeros(Nw,1);
                    for pp = 1:Nw
                        aa = (bwater == pp);
                        waterSize(pp) = sum(aa(:));
                    end
                    [val, ind] = max(waterSize);
                    % assign this cluster 1, the rest (dry areas) is 0
                    bb = (bwater == ind);
                    gridcell(bb) = 1;
                    gridcell(~bb) = 0;
                else
                    % if the entire cell is dry, assign 0 to all pixels
                    gridcell(:) = 0;
                end
                % ------ step 2: check inter-cell connectivity ------
                if intercell == 1
                    fNm = (gridcell(1,:) == 1);
                    fNp = (gridcell(end,:) == 1);
                    fOm = (gridcell(:,1) == 1);
                    fOp = (gridcell(:,end) == 1);
                    if sum(fNm) == 0
                        block.Nm(ii,jj,kk) = 1;
                    end
                    if sum(fNp) == 0
                        block.Np(ii,jj,kk) = 1;
                    end
                    if sum(fOm) == 0
                        block.Om(ii,jj,kk) = 1;
                    end
                    if sum(fOp) == 0
                        block.Op(ii,jj,kk) = 1;
                    end 
                end
                bathBinary((ii-1)*r+1:ii*r,(jj-1)*r+1:jj*r) = gridcell;
            end
        end
        % ------ step 3: check cross-cell connectivity ------
        if crosscell == 1
            for ii = 2:Dim(1)
                for jj = 2:Dim(2)
                    % extract the wet areas of neighboring cells
                    thiscell = bathBinary((ii-1)*r+1:ii*r,(jj-1)*r+1:jj*r);
                    cellNm = bathBinary((ii-2)*r+1:(ii-1)*r,(jj-1)*r+1:jj*r);
                    cellOm = bathBinary((ii-1)*r+1:ii*r,(jj-2)*r+1:(jj-1)*r);
                    % check if the wet pixels on adjacent cell faces are connected
                    % the x-face
                    fCNm = (thiscell(1,:) == 1);
                    fMNm = (cellNm(end,:) == 1);
                    flagN = 0;
                    for pp = 1:r
                        if (fCNm(pp) == 1 && fMNm(pp) == 1)
                            flagN = 1;
                            break;
                        end
                    end
                    if flagN == 0
                        block.Nm(ii,jj,kk) = 1;
                        block.Np(ii-1,jj,kk) = 1;
                    end
                    % the y-face
                    fCOm = (thiscell(:,1) == 1);
                    fMOm = (cellOm(:,end) == 1);
                    flagO = 0;
                    for pp = 1:r
                        if (fCOm(pp) == 1 && fMOm(pp) == 1)
                            flagO = 1;
                            break;
                        end
                    end
                    if flagO == 0
                        block.Om(ii,jj,kk) = 1;
                        block.Op(ii,jj-1,kk) = 1;
                    end


                end
            end
        end
    end

end

