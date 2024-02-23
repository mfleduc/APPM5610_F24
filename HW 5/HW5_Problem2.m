clear all;close all;clc;
rng(45735684)
%% HW 5 Problem 2
e = 10^-6;
A = [[2,e];[e,1]];
evalsNoShift = QRIteration(A,1,10^-8);%No shift
evalsShift = QRIteration(A,1,10^-8,1);%Shift mu=1
evals = eig( A );

diffNoShift = max( abs(evals-sort(evalsNoShift)) );
diffShift =  max( abs(evals-sort(evalsShift)) );

%% HW 5 Problem 3
ns = [6, 12, 18, 24];
Bs = {};
for n=ns
    B = diag( 1:n )+diag( rand(1,n-1),1 ) ;
    B=B+B' ;
    Bs{end+1}=B;
    %Constructing arbitrary tridiagonal matrix with eigenvalues approximately
    %2:2:2n (by Gershgorin Circle Thm)
    %No shift first:
    nIters = 1000;
    tol = 1e-8;
    [evalsB,errCode,dels] = QRIteration(B,nIters,tol);
    %Shifted QR, mu=1
    [evalsBs,errCodes,delss] = QRIteration(B,nIters,tol,1);
    %Shifted QR, mu=-1
    [evalsBn,errCoden,delsn] = QRIteration(B,nIters,tol,-1);
    KB = length(dels);KBs = length(delss);KBn = length(delsn);
    evalsBs = sort(evalsBs);
    evalsB = sort(evalsB);
    figure;
    semilogy( 1:KB, dels,'b','linewidth',2  );
    grid on;
    hold on;
    semilogy( 1:KBs, delss, 'red','linewidth',2 );
    semilogy(1:KBn,delsn, 'g','linewidth',2);
%     semilogy( 1:KB,(evalsB(end-1)/evalsB(end)).^(1:KB),'k.', 'linewidth',2)
    xlabel('Number of iterations')
    ylabel('Error')
    title(sprintf('Shifted and unshifted QR, tridiagonal matrix, n=%d',n))
    legend( 'Unshifted', '\mu=1', '\mu=-1', 'location', 'southwest') 
    savefig( sprintf('convergence_n%d.fig', n ));
    saveas(gcf, sprintf('convergence_n%d.png', n ));
end
