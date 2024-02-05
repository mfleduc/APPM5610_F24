function [lambda, eigvec,errMsg,cnt] = powerMethod(A,tol,nIters,varargin)
%% powerMethod(A,tol,nIters,varargin)
%Implementation of the power method for finding the largest (in absolute value)
%eigenvalue of the matrix A 
%Inputs:
%   A: The matrix to find the larest eigenvalue of
%   tol: Relative error tolerance for convergence
%   nIters: Maximum nmber of iterations
%   varargin{1}: Initial guess of the eigenvector. If not supplied, will
%       start with a random vector. 
%Outputs:
%   lambda: Estimate of the eigenvalue
%   eigvec: Estimate of the eigenvector
%   errMsg: Error message
%           =0: Successful run
%           =1: Not successful. Reached maximum numer of iterations

[n,m] = size(A);
if n ~= m
    error('A must be a square matrix');
end
if isempty(varargin)%Generating random vector
    eigvec = randn(n,1);
else
    eigvec = varargin{1};
end
eigvec = eigvec/norm(eigvec);%Normalize the initial guess
del = tol+1;
cnt=0; 
lambda = 1;
while del>tol && cnt<nIters
    cnt=cnt+1;
    zk = A*eigvec;%Stepping
    qk = zk/norm(zk);%Normalizing
    lk = sum(qk.*(A*qk));%Inner product
    dL = abs(lk-lambda)/abs(lk);%Calculating relative errors
    dV = norm(qk-eigvec);
    del = max(dL,dV);
    lambda = lk;
    eigvec = qk;
end
if cnt<nIters
    errMsg=0;
else
    errMsg=1;
end
end

