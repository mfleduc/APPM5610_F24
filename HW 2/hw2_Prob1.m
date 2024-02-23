%% APPM 5610 HW 2
clear variables;close all;clc;
rng(480958)
%%Part B
tol = 1e-6;
maxIters = 1e5;
ns = 2.^(1:12);

evecBiggest = {};
lambdaBiggest = [];
fprintf('           PART B \n');
figure;
legendStr = {};
for nn = 1:length(ns)
H = hilb(ns(nn));
[lambdaBiggest(nn),evecBiggest{nn}] = powerMethod( H,tol,maxIters );
fprintf('n=%d: The largest eigenvalue in absolute value is approximately %.4f \n',ns(nn), lambdaBiggest(nn));
plot( (0:ns(nn)-1)/(ns(nn)-1), evecBiggest{nn}*sign(evecBiggest{nn}(1)) );hold on;
legendStr{nn} = sprintf('n=%d', ns(nn));
end
grid on;
xlabel('Entries of the eigenvector')
ylabel('Value of entry')
legend(legendStr,'location', 'northeast')
ylim([0 2])
title('Dominant eigenvectors of the nxn Hilbert matrix. n=2^{(1:12)}')
%%Part C
%%Shifted power method
fprintf('           PART C \n');
A = H-lambdaBiggest*eye(16); %Should allow us to get the smallest eigenvalue
[lambdaSmallest,evecSmallest,errMsgSmallest,cntsmall] = powerMethod( A,tol,maxIters );
fprintf('Using the shifted method the smallest eigenvalue in absolute value is approximately %.4f \n', lambdaSmallest+lambdaBiggest);
%%Inverse power method
Hinv = invhilb(16);
[laminv,evecinv,errMsginv,cntinv] = powerMethod( Hinv,tol,maxIters );
fprintf('Using the inverse method the smallest eigenvalue in absolute value is approximately %.4f \n', 1./laminv);
%% Testing idea for part D
e = 0.1;
D = diag([1,-1,0.5]);
[eval1,evec1] = powerMethod(D+e*eye(3),tol,maxIters);
[eval2,evec2] = powerMethod(D-e*eye(3),tol,maxIters);
fprintf('           PART D         \n');
fprintf('Negative eigenvalue is approximately %.4f, eigenvector is (%.4f, %.4f, %.4f)\n',eval2+e,evec2);
fprintf('Positive eigenvalue is approximately %.4f, eigenvector is (%.4f, %.4f, %.4f)\n', eval1-e, evec1);
%% Testing idea for part E
L = diag([1,1,1,1,1,1,1,1,1/2,1/4]);%Standard unit vectors are the eigenvectors
evecGuesses = zeros(10,8);
evalGuesses = [];
cnts=[];
for kk = 1:8
    [evalGuesses(kk), evecGuesses(:,kk),~,cnts(kk)] = powerMethod(L, tol,maxIters);
end
tmp = evecGuesses(1:8,:)\eye(8);
fprintf('           PART E       \n');
fprintf('The calculated eigenvectors are \n');
disp(evecGuesses)
fprintf('Checking if the standard unit vectors are in the span of the results from the power method.\n');
fprintf('Expect the columns to be the standard unit vectors in R^8: evecGuesses*tmp =  \n');
disp(evecGuesses(1:8,:)*tmp);
