clear all; close all;

% Adjustable parameters.
J = 16
aa = 1.3

% Define the test function
xx = (-pi:2*pi/J:pi-2*pi/J)';
bb = sqrt(8/(pi*(1+16*aa^2)));
uu = bb*( (aa+cos(xx)).^2 - (0.5+aa^2) );
uk = fft(uu)
% Exact Hamiltonian values
h2ex = 1 + 3/(1+16*aa^2);
h3ex = 8*sqrt(2)*pi*aa^2 * (16*pi*aa^2 + pi)^(-3/2);
% Compute the Hamiltonian with the routine
[ham,ham2,ham3] = hamiltonian(uk,1.,1.,1.);

ham2
h2ex
ham3
h3ex

h2relerr = abs(h2ex-ham2)/abs(h2ex)
h3relerr = abs(h3ex-ham3)/abs(h3ex)