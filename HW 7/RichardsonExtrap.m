function out = RichardsonExtrap(calcedVals,hRatio,tol)
%% out = RichardsonExtrap(in)
%Implementation of Ricahrdson Extrapolation
%Inputs:
%   calcedVals: Vector of values computed for various h. The first entry
%   should be calculated with the largest step size, second with second
%   largest, etc.
%   hRatio: Ratio of the step size used for each entry in calcedVals to the
%   step size used to calculate the first value. 
%   tol: relative error tolerance. Stops when there is a calculated value
%   within tol OR you run out of room to iterate.
N = length(calcedVals);
R = zeros(N);R(:,1) = calcedVals;

for jj = 1:N-1
    for ii = 1:N
        keyboard
    end
end

out = calcedVals(end);
end