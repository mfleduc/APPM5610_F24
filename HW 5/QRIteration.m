function [evals,errCode,delVec] = QRIteration(A,nIters,tol,varargin)
%% [evals,errCode] = QRIteration(A,nIters,tol,varargin)
%   Implementation of QR iteration for calculating eigenvalues. Can be
%   shifted or not. 
%Inputs:
%   A: Square matrix whose eigenvalues are desired.
%   nIters: Maximum number of iterations allowed
%   tol: Tolerance desired for convergence. Convergence is achieved when
%   off(A)/||diag(A)|| < tol, where off(A) = ||A-diag(A)||_{fro}
%Outputs:
%   evals: Vector of eigenvalues
%   errCode: Error code. 0 if everything went well, 1 if the code failed to
%   converge in nIters.
%   delVec: For testing purposes, tracking the convergence of the algorithm

del = tol+1;
delVec = zeros( 1,nIters );
n = size(A,1);
if isempty(varargin)
    mu=0;%No shifting, but we still need a value of \mu to correct the evals later
else
    
    mu = varargin{1};
    A = A+mu*eye(n);%Shift the iteration 
end
%Iterate until off(A)<tol
cnt = 0;
while cnt<nIters && del>tol
    [Q,R] = qr(A);
    A=R*Q;
    del = norm(A-diag(diag(A)),'fro')/norm(diag(A));%Frobenius norm of the off diagonal terms
    cnt=cnt+1;
    delVec(cnt)=del;
end
evals = diag(A)-mu;
delVec = delVec(1:cnt);
if cnt == nIters
    errCode=1;
else
    errCode=0;
end

end

