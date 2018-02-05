function [subA] = ComputeSubgridFaceArea(subA, depth, dim, dx, r, kk)

% Given a pre-defined free surface surf(kk), this function computes the
% subgrid face areas for each coarse grid cell.

Dim = dim / r;

for ii = 1:Dim(1)
    for jj = 1:Dim(2)
        i1 = (ii-1)*r+1;
        i2 = ii*r;
        j1 = (jj-1)*r+1;
        j2 = jj*r;
        % compute subA.Np and subA.Nm
        if i2 ~= dim(1)
            depthNp = depth(i2,j1:j2);
            depthNm = depth(i2+1,j1:j2);
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
            depthOp = depth(i1:i2,j2);
            depthOm = depth(i1:i2,j2+1);
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



end