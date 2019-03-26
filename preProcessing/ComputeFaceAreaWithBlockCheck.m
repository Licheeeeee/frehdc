function [subA, block] = ComputeFaceAreaWithBlockCheck(subA, bathCenter, surf, r, dx)
%% ZHiLi 20190129

dim = size(bathCenter);
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
    map.Np = zeros(Dim(1),Dim(2),r/2,r);
    map.Nm = zeros(Dim(1),Dim(2),r/2,r);
    map.Op = zeros(Dim(1),Dim(2),r,r/2);
    map.Om = zeros(Dim(1),Dim(2),r,r/2);
    fprintf('Surface elevation = %f ...\n',surf(kk));        
    h = surf(kk);
    depth = max(h - bathCenter, 0);
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
                % remove small water area, ZhiLi20180214
                if (sum(bb(:)) < 0.01*r*r)
                    gridcell(:) = 0;
                else
                    gridcell(bb) = 1;
                    gridcell(~bb) = 0;
                end
            else
                % if the entire cell is dry, assign 0 to all pixels
                gridcell(:) = 0;
            end
            map.Np(ii,jj,:,:) = gridcell(r/2+1:end,:);
            map.Nm(ii,jj,:,:) = gridcell(1:r/2,:);
            map.Op(ii,jj,:,:) = gridcell(:,r/2+1:end);
            map.Om(ii,jj,:,:) = gridcell(:,1:r/2);
            % ------ step 2: check inter-cell connectivity ------
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
            bathBinary((ii-1)*r+1:ii*r,(jj-1)*r+1:jj*r) = gridcell;
        end
    end
    % ------ step 3: check cross-cell connectivity ------
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
    % -------- step 4: compute face areas ------
    for ii = 1:Dim(1)
        for jj = 1:Dim(2)
            i1 = (ii-1)*r+1;
            i2 = ii*r;
            j1 = (jj-1)*r+1;
            j2 = jj*r;
            % compute subA.Np and subA.Nm
            if i2 ~= dim(1)
                % Np face of the current cell
                depthNp = depth(i2,j1:j2);
                mapNp = squeeze(map.Np(ii,jj,end,:));
                aa = mapNp == 0;
                depthNp(aa) = 0;
                % Nm face of the iP cell
                depthNm = depth(i2+1,j1:j2);
                mapNm = squeeze(map.Nm(ii+1,jj,1,:));
                aa = mapNm == 0;
                depthNm(aa) = 0;
                % compute mutually-wet area
                subA.Np(ii,jj,kk) = sum(min(depthNp, depthNm)) * dx;
                subA.Nm(ii+1,jj,kk) = sum(min(depthNp, depthNm)) * dx;
            else
                depthNp = depth(i2,j1:j2);
                subA.Np(ii,jj,kk) = sum(depthNp) * dx;
                depthNm = depth(1,j1:j2);
                subA.Nm(1,jj,kk) = sum(depthNm) * dx;
            end
            % compute subA.Op and subA.Om
            if j2 ~= dim(2)
                % Op face of the current cell
                depthOp = depth(i1:i2,j2);
                mapOp = squeeze(map.Op(ii,jj,:,end));
                aa = mapOp == 0;
                depthOp(aa) = 0;
                % Om face of the iP cell
                depthOm = depth(i1:i2,j2+1);
                mapOm = squeeze(map.Om(ii,jj+1,:,1));
                aa = mapOm == 0;
                depthOm(aa) = 0;
                % compute mutually-wet areas
                subA.Op(ii,jj,kk) = sum(min(depthOp, depthOm)) * dx;
                subA.Om(ii,jj+1,kk) = sum(min(depthOp, depthOm)) * dx;
            else
                depthOp = depth(i1:i2,j2);
                subA.Op(ii,jj,kk) = sum(depthOp) * dx;
                depthOm = depth(i1:i2,1);
                subA.Om(ii,1,kk) = sum(depthOm) * dx;
            end
        end
    end
    
    
    for ii = 1:Dim(1)
        for jj = 1:Dim(2)
            i1 = (ii-1)*r+1;
            i2 = ii*r;
            j1 = (jj-1)*r+1;
            j2 = jj*r;

            % compute min and max areas
            if i2 ~= dim(1)
                depthN = depth(i2-floor(r/2)+1:i2+round(r/2),j1:j2);
                depthNex = depth(i2-floor(r*3/4)+1:i2+round(r*3/4),j1:j2);
                mapN = squeeze(cat(3,map.Np(ii,jj,:,:),map.Nm(ii+1,jj,:,:)));
            else
                depthN = depth(i2-floor(r/2)+1:i2,j1:j2);
                mapN = squeeze(map.Np(ii,jj,:,:));
            end
%             aa = mapN == 0;
%             depthN(aa) = 0;
            sumN = sum(depthN,2);
            sumNex = sum(depthNex,2);
            bb = sumN == 0;
            sumN(bb) = [];
            if ~isempty(sumN)
                Nmin = min(sumN);
                Nminex = min(sumNex);
                Nmax = max(sumN);
                Nmed = median(sumN);
                Nmid = subA.Np(ii,jj,kk);
                subA.Nmin(ii,jj,kk) = Nmin;
                
                % find location of min area
                cc = sum(depthNex,2) == Nminex;
                % compute distance from minA to cell face
                ind = find(cc == 1);
                dist = min([abs(ind(1)-r/2) abs(ind(end)-r/2)]);
                if i2 ~= dim(1) && bb(1) == 0 && bb(end) == 0 && ...
                        subA.Np(ii+1,jj,kk) > 0 && subA.Nm(ii,jj,kk) > 0 &&...
                        length(sumN) > 0.5*r && cc(1) == 0 && cc(end) == 0
                    % reduce area for internal contraction

                    if length(sumN) > 3

                        if (abs(Nmax-Nmed) > 2*abs(Nmed-Nmin))
                            subA.Nmax(ii,jj,kk) = Nmin;% + 2*dist*(Nmid-Nmin)/r;
                        elseif (abs(Nmed-Nmin) > 2*abs(Nmax-Nmed))
                            subA.Nmax(ii,jj,kk) = Nmin;% + 2*dist*(Nmid-Nmin)/r;
                        else
                            subA.Nmax(ii,jj,kk) = Nmid;
                        end

                    else
                        subA.Nmax(ii,jj,kk) = Nmid;
                    end
                else
                    subA.Nmax(ii,jj,kk) = Nmid;
                end
            end
            
            
%             if Nmin ~= Nmax
%                 ismin = sumN == Nmin;
%                 % if min area not on boundary, use min
%                 if ismin(1) == 0 && ismin(end) == 0
%                     subA.Nmin(ii,jj,kk) = Nmin;
%                     subA.Nmax(ii,jj,kk) = Nmax;
%                 else
%                     % min area on boundary, mid area equals max or min
%                     if Nmid == Nmin || Nmid == Nmax
%                         subA.Nmin(ii,jj,kk) = Nmin;
%                         subA.Nmax(ii,jj,kk) = Nmax;
%                     else
%                         subA.Nmin(ii,jj,kk) = Nmid;
%                         subA.Nmax(ii,jj,kk) = Nmid;
%                     end
%                 end
%             else
%                 subA.Nmin(ii,jj,kk) = Nmin;
%                 subA.Nmax(ii,jj,kk) = Nmax;
%             end
            
            
            
            if j2 ~= dim(2)
                depthO = depth(i1:i2,j2-floor(r/2)+1:j2+round(r/2));
                depthOex = depth(i1:i2,j2-floor(3*r/4)+1:j2+round(3*r/4));
                mapO = squeeze(cat(4,map.Op(ii,jj,:,:),map.Om(ii,jj+1,:,:)));
            else
                depthO = depth(i1:i2,j2-floor(r/2)+1:j2);
                mapO = squeeze(map.Op(ii,jj,:,:));
            end
            aa = mapO == 0;
            depthO(aa) = 0;
            sumO = sum(depthO,1);
            sumOex = sum(depthOex,1);
            bb = sumO == 0;
            sumO(bb) = [];
%             if ii == 2 && jj == 84
%                 disp('www');
%             end
            if ~isempty(sumO)
                Omin = min(sumO);
                Ominex = min(sumOex);
                Omax = max(sumO);
                Omed = median(sumO);
                Omid = subA.Op(ii,jj,kk);
                subA.Omin(ii,jj,kk) = Omin;
                cc = sum(depthOex,1) == Ominex;
                % compute distance from minA to cell face
                ind = find(cc == 1);
                dist = min([abs(ind(1)-r/2) abs(ind(end)-r/2)]);
                % only reduce area if this cell is fully connected
                if j2 ~= dim(2) && bb(1) == 0 && bb(end) == 0 && ...
                        subA.Op(ii,jj+1,kk) > 0 && subA.Om(ii,jj,kk) > 0 &&...
                        length(sumO) > 0.5*r && cc(1) == 0 && cc(end) == 0
                        
                    if length(sumO) > 3
                        if (abs(Omax-Omed) > 2*abs(Omed-Omin))
                            subA.Omax(ii,jj,kk) = Omin;% + 2*dist*(Omid-Omin)/r;
                        elseif (abs(Omed-Omin) > 2*abs(Omax-Omed))
                            subA.Omax(ii,jj,kk) = Omin;% + 2*dist*(Omid-Omin)/r;
                        else
                            subA.Omax(ii,jj,kk) = Omid;
                        end
                    else
                        subA.Omax(ii,jj,kk) = Omid;
                    end
                else
                    subA.Omax(ii,jj,kk) = Omid;
                end
            end
            
%             if Omin ~= Omax
%                 ismin = sumO == Omin;
%                 % if min area not on boundary, use min
%                 if ismin(1) == 0 && ismin(end) == 0
%                     subA.Omin(ii,jj,kk) = Omin;
%                     subA.Omax(ii,jj,kk) = Omax;
%                 else
%                     % min area on boundary, mid area equals max or min
%                     if Omid == Omin || Omid == Omax
%                         subA.Omin(ii,jj,kk) = Omin;
%                         subA.Omax(ii,jj,kk) = Omax;
%                     else
%                         subA.Omin(ii,jj,kk) = Omid;
%                         subA.Omax(ii,jj,kk) = Omid;
%                     end
%                 end
%             else
%                 subA.Omin(ii,jj,kk) = Omin;
%                 subA.Omax(ii,jj,kk) = Omax;
%             end
        end
    end
end