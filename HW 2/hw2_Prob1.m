%% APPM 5610 HW 2
clear variables;close all;clc;
rng(480958)
%Part B
tol = 1e-8;
maxIters = 1e6;
H = hilb(16);
[lambdaBiggest,evecBiggest,errMsgBiggest,cntbig] = powerMethod( H,tol,maxIters );
fprintf('The largest eigenvalue in absolute value is approximately %.4f \n', lambdaBiggest);
%Part C
%Shifted power method
A = H-lambdaBiggest*eye(16); %Should allow us to get the smallest eigenvalue
[lambdaSmallest,evecSmallest,errMsgSmallest,cntsmall] = powerMethod( A,tol,maxIters );
fprintf('Using the shifted method the smallest eigenvalue in absolute value is approximately %.4f \n', lambdaSmallest+lambdaBiggest);
%Inverse power method
Hinv = invhilb(16);
[laminv,evecinv,errMsginv,cntinv] = powerMethod( Hinv,tol,maxIters );
fprintf('Using the inverse method the smallest eigenvalue in absolute value is approximately %.4f \n', 1./laminv);