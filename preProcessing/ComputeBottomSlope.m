%% Compute bottom slope due to macro structures, ZhiLi20180417
function [cv] = ComputeBottomSlope(bathCenter, surf, r)
    % ---------- User Settings ----------
    % bathCenter = fine grid data
    % surf = all possible free surface elevations
    % dx = the grid size hierachy
    %      the first element is the input fine grid size
    %      the last element is the output coarse grid size
    
    dx = 1;
    % the input data size
    dim = size(bathCenter);
    % the output data size
    Dim = dim / r;    
    N = length(surf);
    % initialize output structure
    cv.CvY = zeros([Dim N]);
    cv.CvX = zeros([Dim N]);
    % ------------ Computation ------------
    % loop over all possible surface elevations
    for kk = 1:N
        fprintf('Bottom Slope ---> Surface elevation = %f ...\n',surf(kk));        
        % ------ compute the bottom curvature ------
        % added by ZhiLi 20180409
        % in y direction
        for ii = 1:Dim(1)
            for jj = 1:Dim(2)-1
                % extract one grid cell centered at plus face
                gridcell = bathCenter((ii-1)*r+1:ii*r,(jj-1)*r+1:(jj+1)*r);
                depth = max(surf(kk) - gridcell, 0);                
                CvY = 0;
                % compute slope within each coarse cell
                % uphill in positive direction is positive slope -->
                % increase drag coefficient
                for qq = 1:r
                    for ll = 1:2*r-1
                        if depth(qq,ll) > 0 && depth(qq,ll+1) > 0
                            thisCv = (depth(qq,ll+1) - depth(qq,ll)) ./ dx ./ (ll-(2*r+1)/2);
                            CvY = CvY + thisCv / (depth(qq,ll)^(1/3));
                        end
                    end
                end
                cv.CvY(ii,jj,kk) = CvY;              
            end
        end
        % in x direction
        for ii = 1:Dim(1)-1
            for jj = 1:Dim(2)
                % extract one grid cell centered at plus face
                gridcell = bathCenter((ii-1)*r+1:(ii+1)*r,(jj-1)*r+1:jj*r);
                depth = max(surf(kk) - gridcell, 0);                
                CvX = 0;
                % compute curvature within each coarse cell
                for ll = 1:r
                    for qq = 2:2*r-1
                        if depth(qq,ll) > 0 && depth(qq+1,ll) > 0
                            thisCv = (depth(qq+1,ll) - depth(qq,ll)) ./ dx ./ (qq-(2*r+1)/2);
                            CvX = CvX + thisCv / (depth(qq,ll)^(1/3));
                        end
                    end
                end
                cv.CvX(ii,jj,kk) = CvX; 
            end
        end
    end

end

