%% Problem 2
clear variables
close all
t = 0:0.0001:(2*pi);
zAB2 = 2*(exp(2*1j*t)-exp(1j*t))./(3*exp(1j*t)-1);
zAB3 = 12*(exp(3*1j*t)-exp(2j*t))./(23*exp(2j*t)-16*exp(1j*t)+5);
zAM2 = 12*(exp(3*1j*t)-exp(2j*t))./(5*exp(2j*t)+8*exp(1j*t)-1) ;

figure; 
plot(real(zAB2),imag(zAB2), 'b','linewidth',2)
hold on;grid on;
plot(real(zAB3),imag(zAB3), 'r','linewidth',2)
plot(real(zAM2),imag(zAM2), 'k','linewidth',2)
legend('2nd order Adams-Bashforth','3rd order Adams-Bashforth', '2nd order Adams-Moulton')
title('Regions of absolute stability for the given schemes')
xlabel('Re(z)')
ylabel('Im(z)')

%% Problem 3
zBD2 = (exp(2*1j*t)-4/3*exp(1j*t)+1/3)./(2/3*exp(2*1j*t));
zBD3 = (exp(3*1j*t)-18/11*exp(2*1j*t)+9/11*exp(1j*t)-2/11)./(6/11*exp(3*1j*t));

figure; 
plot(real(zBD2),imag(zBD2), 'b','linewidth',2)
hold on;grid on;
plot(real(zBD3),imag(zBD3), 'r','linewidth',2)

legend('2 step back-diff','3 step back-diff')
title('Regions of absolute stability for the given schemes')
xlabel('Re(z)')
ylabel('Im(z)')