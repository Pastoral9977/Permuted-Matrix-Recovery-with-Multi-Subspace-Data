%--------------------------------------------------------------------------
% This is the function to call the sparse optimization program, to call the 
% spectral clustering algorithm and to compute the clustering error.
% r = projection dimension, if r = 0, then no projection
% affine = use the affine constraint if true
% s = clustering ground-truth
% missrate = clustering error
% CMat = coefficient matrix obtained by SSC
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function [grps] = SSC_with_numgroups_only(X,n)


    rho = 1;
    outlier = false;
    alpha = 20;
    affine = false;
    r = 0;

    Xp = DataProjection(X,r);

    if (~outlier)
        CMat = admmLasso_mat_func(Xp,affine,alpha);
        C = CMat;
    else
        CMat = admmOutlier_mat_func(Xp,affine,alpha);
        N = size(Xp,2);
        C = CMat(1:N,:);
    end

    CKSym = BuildAdjacency(thrC(C,rho));
    grps = SpectralClustering(CKSym,n);

end