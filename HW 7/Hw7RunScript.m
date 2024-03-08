%% LeDuc APPM 5610 HW 7 run script
h = 10^-2;
factor = 0.5.^(0:4);
rhs = @Hw7Rhs;
y0 = [1/2;0];
tint = [0,3*pi];
out = TrapezoidRule(rhs, y0,h,tint);
y = out(1,:).*t;
vals = zeros(length(factor),1);
for kk = 1:lenth(factor)
    out = TrapezoidRule(rhs, y0,h,tint);
    vals(kk) = out(1,end)*3*pi;
end
reVal = RichardsonExtrap(vals,factor, 1e-10) ;
