%% LeDuc APPM 5610 HW 7 run script
%Problem 1
h = 10^-3;
tol=1e-10;
factor = 0.5.^(1:5);
rhs = @Hw7Rhs;
y0 = [1/2;0];
tint = [0,3*pi];
vals = zeros(length(factor),1);
% for kk = 1:length(factor)
%     out = TrapezoidRule(rhs, y0,h*factor(kk),tint);
%     vals(kk) = out(1,end)*3*pi;
% end
%Richardson extrapolation
prefactor = 2;
N = length(factor)+1;
R = zeros(N);
trueVal = besselj(1,3*pi);
out = TrapezoidRule(rhs, y0,h,tint);
R(1,1) = out(1,end)*3*pi;
for ii = 1:N-1
    out = TrapezoidRule(rhs, y0,h*factor(ii),tint);
    R(ii+1,1) = out(1,end)*3*pi;
    for jj = 1:ii
        R(ii+1,jj+1) = ((2^(jj+1))*R(ii + 1, jj) - R(ii, jj))/(2^(jj+1) - 1);
    end
    del = abs(R(ii+1,ii+1)-trueVal)/abs(trueVal) ;
    if del<tol
       break 
    end
end