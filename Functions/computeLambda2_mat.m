%--------------------------------------------------------------------------
% This function takes a DxN matrix of N data points in a D-dimensional 
% space and returns the regularization constant of the L1 norm
% Y: DxN data matrix
% lambda: regularization parameter for lambda*||C||_1 + 0.5 ||Y-YC||_F^2
%--------------------------------------------------------------------------
% Copyright @ Ehsan Elhamifar, 2012
%--------------------------------------------------------------------------

function [lambda,value] = computeLambda2_mat(Y,P)

if (nargin < 2)
    P = Y;
end

N = size(Y,2);
T = P' * Y;
T(1:N,:) = T(1:N,:) - diag(diag(T(1:N,:)));
T = abs(T);
value = max(T,[],1);
lambda = min(value)*0.8;
lambda = 1/lambda;

end



