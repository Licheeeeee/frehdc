%% Compute drag correction due to macro structures, ZhiLi20180327
function [cd] = ComputeDragCorrection(bathCenter, surf, r)
    % ---------- User Settings ----------
    % bathCenter = fine grid data
    % surf = all possible free surface elevations
    % dx = the grid size hierachy
    %      the first element is the input fine grid size
    %      the last element is the output coarse grid size
    
    % the input data size
    dim = size(bathCenter);
    % the output data size
    Dim = dim / r;    
    N = length(surf);
    % initialize output structure
    cd.effCdY = zeros([Dim N]);
    cd.effCdX = zeros([Dim N]);
    % ------------ Computation ------------
    % loop over all possible surface elevations
    for kk = 1:N
        fprintf('Sidewall Drag ---> Surface elevation = %f ...\n',surf(kk));        
        aa = (surf(kk) > bathCenter);
        bathBinary = zeros(dim);
        % wet area is 1, dry area is 0
        bathBinary(aa) = 1;
        bathBinary(~aa) = 0;
        % ------ compute the effective drag ------
        % added by ZhiLi 20180327
        % in y direction
        for ii = 1:Dim(1)
            for jj = 1:Dim(2)-1
                % extract one grid cell centered at plus face
                gridcell = bathBinary((ii-1)*r+1:ii*r,(jj-1)*r+1:(jj+1)*r);
                aa = gridcell == 1;
                % fine the object boundaries for partially wet cells
                if sum(aa(:)) < 2*r^2 && sum(aa(:)) > 0
                    % set water=0 and land=1, find objects in the interior
                    % of the cell
                    gridcell(aa) = 0;
                    gridcell(~aa) = 1;
                    % sketch the outline of land
                    gridcell2 = zeros(size(gridcell));
                    for gii = 1:size(gridcell,1)
                        for gjj = 1:size(gridcell,2)
                            if gridcell(gii,gjj) == 1
                                gridcell2(gii,gjj) = 1;
                                if gii > 1
                                    gridcell2(gii-1,gjj) = 1;
                                end
                                if gii < size(gridcell,1)
                                    gridcell2(gii+1,gjj) = 1;
                                end
                                if gjj > 1
                                    gridcell2(gii,gjj-1) = 1;
                                end
                                if gjj < size(gridcell,2)
                                    gridcell2(gii,gjj+1) = 1;
                                end
                            end
                        end
                    end
                    [Bind, Lmap, Nb] = bwboundaries(gridcell2);
                    % loop over all boundaries
                    if Nb > 0
%                         if jj == 11
%                             [ii jj]
%                             disp('pause');
%                         end
                        for ll = 1:Nb
                            thisB = Bind{ll,1};
                            thisDr = 0;
                            % thisBin is the object boundary cells in the
                            % interior of the coarse grid cell
                            for qq = 2:size(thisB,1)
                                if thisB(qq,1) ~= thisB(qq-1,1) && thisB(qq,1) ~= 1 && thisB(qq,1) ~= size(gridcell2,1)
                                    thisDr = thisDr + 1;% / abs(thisB(qq,2) - (2*r+1)/2);
                                end
                            end
                            % if the object has horizontal span, the drag
                            % is reduced
                            stdB = std(thisB(:,2));
                            if stdB > 0
                                thisDr = thisDr / (stdB+1);
                            end
                            cd.effCdY(ii,jj,kk) = cd.effCdY(ii,jj,kk) + thisDr;
                        end
                    end
                end
            end
        end
        % in x direction
        for ii = 1:Dim(1)-1
            for jj = 1:Dim(2)
                % extract one grid cell centered at plus face
                gridcell = bathBinary((ii-1)*r+1:(ii+1)*r,(jj-1)*r+1:jj*r);
                aa = gridcell == 1;
                % fine the object boundaries
                if sum(aa(:)) < 2*r^2 && sum(aa(:)) > 0
                    % set water=0 and land=1, find objects in the interior
                    % of the cell
                    gridcell(aa) = 0;
                    gridcell(~aa) = 1;
                    % sketch the outline of land
                    gridcell2 = zeros(size(gridcell));
                    for gii = 1:size(gridcell,1)
                        for gjj = 1:size(gridcell,2)
                            if gridcell(gii,gjj) == 1
                                gridcell2(gii,gjj) = 1;
                                if gii > 1
                                    gridcell2(gii-1,gjj) = 1;
                                end
                                if gii < size(gridcell,1)
                                    gridcell2(gii+1,gjj) = 1;
                                end
                                if gjj > 1
                                    gridcell2(gii,gjj-1) = 1;
                                end
                                if gjj < size(gridcell,2)
                                    gridcell2(gii,gjj+1) = 1;
                                end
                            end
                        end
                    end
                    [Bind, Lmap, Nb] = bwboundaries(gridcell2);
                    % loop over all boundaries
                    if Nb > 0
                        for ll = 1:Nb
                            thisB = Bind{ll,1};
                            thisDr = 0;
                            % thisBin is the object boundary cells in the
                            % interior of the coarse grid cell
                            for qq = 2:size(thisB,1)
                                if thisB(qq,2) ~= thisB(qq-1,2) && thisB(qq,2) ~= 1 && thisB(qq,2) ~= size(gridcell2,2)
                                    thisDr = thisDr + 1;% / abs(thisB(qq,1) - (2*r+1)/2);
                                end
                            end
                            % if the object has horizontal span, the drag
                            % is reduced
                            stdB = std(thisB(:,1));
                            if stdB > 0
                                thisDr = thisDr / (stdB+1);
                            end
                            cd.effCdX(ii,jj,kk) = cd.effCdX(ii,jj,kk) + thisDr;
                        end
                    end
                end
            end
        end
    end

end

