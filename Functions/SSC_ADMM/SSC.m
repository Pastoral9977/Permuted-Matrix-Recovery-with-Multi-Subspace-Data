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

function [missrate,ggrps] = SSC(X,r,affine,alpha,outlier,rho,s,broadcast)

    if nargin < 8
        broadcast = true;
    end
    
    if (nargin < 6)
        rho = 1;
    end
    if (nargin < 5)
        outlier = false;
    end
    if (nargin < 4)
        alpha = 20;
    end
    if (nargin < 3)
        affine = false;
    end
    if (nargin < 2)
        r = 0;
    end

    n = length(unique(s));
    Xp = DataProjection(X,r);

    if broadcast
        fprintf('\t\tConstructing adjacency matrix...')
        tic
    end
    
    if (~outlier)
        CMat = admmLasso_mat_func(Xp,affine,alpha);
        C = CMat;
    else
        CMat = admmOutlier_mat_func(Xp,affine,alpha);
        N = size(Xp,2);
        C = CMat(1:N,:);
    end
    
    if broadcast
        fprintf('Done, costing %.3fs\n', toc)

        fprintf('\t\tBuildAdjacency...')
        tic
    end
    
    CKSym = BuildAdjacency(thrC(C,rho));
    if broadcast
        fprintf('Done, costing %.3fs\n', toc)

        fprintf('\t\tSpectralClustering...')
        tic
    end
    grps = SpectralClustering(CKSym,n);
    
    if broadcast
        fprintf('Done, costing %.3fs\n', toc)
    end
    ggrps = bestMap(s,grps);
    missrate = sum(s(:) ~= ggrps(:)) / length(s);
end

