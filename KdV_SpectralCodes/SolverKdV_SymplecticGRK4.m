% This script solves KdV equation in a periodic domain
% u_t=-6uu_x-u_xxx, -L<=x<=L
% using pseudo-spectral & symplectic Generalized Runge-Kutta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set simulation parameters
J = 256;      % Number of points in each direction, J devidable by 4
dt = 1E-3;    % initial time step size
Nt = 10e3;    % Number of time steps
ulim = 1.5E5; % if any u > ulim, simulation stops
Ls = 20;       % half domain size
A=[1/4,1/4-sqrt(3)/6;1/4+sqrt(3)/6,1/4];

% Put useful stuff into a struct
params = struct('Ls',Ls, 'J',J,'dt',dt, 'Nt',Nt,'A',A);
x=(-Ls:2*Ls/J:Ls-2*Ls/J)';

% Set up implicit dispersion
k = (pi/Ls)*[0:J/2 -J/2+1:-1]'; % wavenumbers
Lk = 1i*k.^3;
Dk = 1i*k;

% Initialize
t = 0;
% rng(15); u0 = randn(J,1); 
u0 = 2*sech(x).^2;
uk0 = fft(u0); uk0(J/2+1)=0;
uk = uk0;
ukm1 = uk0;  % at time n-1
ukm2 = uk0;  % at time n-2
% interpolatin function for the initial iteration input
interp = @(u0,um1,um2,w) (.5*w*(w+3)+1).*u0 - w*(w+2).*um1 +.5*w*(w+1).*um2;

figure(10)
plot(x,u0); ylim([-.5 2.5]);
title('solution of KdV equation, symplectic')

% Diagnostics
countDiag = 10;  % Compute diagnostics every countDiag steps
T = zeros(1,Nt/countDiag);
ut = zeros(J,Nt/countDiag);
mass = zeros(1,Nt/countDiag);
energy = zeros(1,Nt/countDiag);
hamiltonian = zeros(1,Nt/countDiag);


% Main loop 
% tic;
for ii=1:Nt
    if mod(ii,countDiag)==0
        if any(isnan(uk(:))), break, end
        T(ii/countDiag) = t;
        ut(:,ii/countDiag) = u;
        
        mass(ii/countDiag) = sum(u)*2*Ls/J;
        energy(ii/countDiag) = sum(abs(uk).^2)*2*Ls/J^2;
        hamiltonian(ii/countDiag) = (real(sum(u2k_dealiasing(uk,params).*conj(uk))) ...
                                     -1/2*sum(abs(Dk.*uk).^2)) *2*Ls/J^2;

%         toc
        display(['iteration i = ', num2str(ii),'; E = ',num2str(energy(ii/countDiag)),...
                 ', H = ',num2str(hamiltonian(ii/countDiag))]);
        figure(10)
        plot(x,u); ylim([-.5 2.5]);
        title('solution of KdV equation, symplectic')
    end
    
    
    % set up initial interation state
    uk10 = interp(uk,ukm1,ukm2,A(1,1)+A(1,2));
    uk20 = interp(uk,ukm1,ukm2,A(2,1)+A(2,2));
    
    % iterations
    [uk1,uk2,ik] = iterativeKdV_GRK4(uk10,uk20, uk,params);
    if ik==1
        disp('Error! quit...');
        break;
    end

    f1 = Lk.*uk1-3*Dk.*u2k_dealiasing(uk1,params);
    f2 = Lk.*uk2-3*Dk.*u2k_dealiasing(uk2,params);
    ukm2 = ukm1;
    ukm1 = uk;
    % Successful step, proceed to evaluation
    uk = uk+.5*dt*(f1+f2);
    t = t+dt;
    u = real(ifft(uk));
    uk = fft(u); uk0(J/2+1)=0;

    if any(abs(u(:))>ulim),break,end
end



% plot results
figure(1)
plot(x,u,'LineWidth',2); hold on;
plot(x,2*sech(mod(x-4*t+Ls,2*Ls)-Ls).^2,'k--');
ylim([-.5 2.5]);
legend('numerics','theory')
title('solution of KdV equation, symplectic');

figure(11)
plot(T,mass); hold on;
title('time-series of the total mass');
figure(12)
plot(T,energy); hold on;
title('time-series of the total energy');
figure(13)
plot(T,hamiltonian); hold on;
title('time-series of the Hamiltonian');

% figure
% contour(x,T,ut',50); colorbar;
% xlabel('x'); ylabel('t');
% title('time-series of the KdV solution')