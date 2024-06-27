%This function aims to divide one picture into four parts for the sake of
%selection

function [X11,X12,X21,X22] = halve(X_tilde,hh,ww)
s = fix(hh/2); t = fix(ww/2); n =size(X_tilde,2);
X11 = zeros(s*t,n);X12 = zeros(s*(ww-t),n);
X21 = zeros((hh-s)*t,n);X22 = zeros((hh-s)*(ww-t),n);

for i = 1:size(X_tilde,2)
    a = X_tilde(:,i);
    X = reshape(a,hh,ww);
    X1 = X(1:s,1:t); 
    X2 = X(1:s,t+1:end); 
    X3 = X(s+1:end,1:t); 
    X4 = X(s+1:end,t+1:end); 
    X11(:,i) = X1(:);
    X12(:,i) = X2(:);
    X21(:,i) = X3(:);
    X22(:,i) = X4(:);
end
end
    
    
    
    