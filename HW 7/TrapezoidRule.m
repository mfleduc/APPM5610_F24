function out = TrapezoidRule(fty,y0, h, tint,varargin)
%% out = trapezoidRule(in)
%Implementation of the trapezoidal rule 
%Inputs:
%   fty: The derivative of y: Right hand side in y'=f(t,y)
%   y0: Initial condition
%   h: Step size
%   tint: Time interval [ t0,tf]

opts = struct();opts.Display='off';%Stopping the code from displaying messages
t = tint(1);
y = y0;
while t<tint(2)
     thisFty = fty(t,y(:,end));
     fn = @(x)(y(:,end)+0.5*h*(thisFty+fty(t+h,x)) - x);
     guess = y(:,end)+h*thisFty;
     yn = fsolve(fn,guess,opts);%Solve the equation that comes from the trapezoidal rule
     t=t+h;
     y(:,end+1) = yn;
end
out = y;
end