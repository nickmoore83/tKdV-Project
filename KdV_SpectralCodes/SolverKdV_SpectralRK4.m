% This script solves KdV equation in a periodic domain
% u_t=-6uu_x-u_xxx, -L<=x<=L
% using pseudo-spectral & explicit and implicit RK4 scheme (non-symplectic)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set simulation parameters
J = 256;      % Number of points in each direction, J devidable by 4
dt = 1E-3;    % initial time step size
Nt = 10e3;    % Number of time steps
ulim = 1.5E5; % if any u > ulim, simulation stops
Ls = 20;       % half domain size

% Put useful stuff into a struct
params = struct('Ls',Ls, 'J',J,'dt',dt, 'Nt',Nt);
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

figure(10)
plot(x,u0); ylim([-.5 2.5]);
title('solution of KdV equation, non-symplectic')

% Diagnostics
countDiag = 10;  % Compute diagnostics every countDiag steps
T = zeros(1,Nt/countDiag);
ut = zeros(J,Nt/countDiag);
mass = zeros(1,Nt/countDiag);
energy = zeros(1,Nt/countDiag);
hamiltonian = zeros(1,Nt/countDiag);

% adaptive stepping stuff:
tol= 1E-1;
r0 = .8*tol;

% Main loop 
% tic;
for ii=1:Nt
    if mod(ii,countDiag)==0
        if any(isnan(uk(:))), break, end
        T(ii/countDiag)=t;
        ut(:,ii/countDiag)=u;
        
        mass(ii/countDiag)=sum(u)*2*Ls/J;
        energy(ii/countDiag) = sum(abs(uk).^2)*2*Ls/J^2; 
%         hamiltonian(ii/countDiag) = sum(u.^3)*2*Ls/J-1/3*sum(abs(Dk.*uk).^2)*2*Ls/J^2;
        hamiltonian(ii/countDiag) = (real(sum(u2k_dealiasing(uk,params).*conj(uk))) ...
                                     -1/2*sum(abs(Dk.*uk).^2)) *2*Ls/J^2;

%         toc
%         display(['iteration i = ', num2str(ii), '; time step dt = ',num2str(dt)]);
        display(['iteration i = ', num2str(ii),'; E = ',num2str(energy(ii/countDiag)),...
                 ', H = ',num2str(hamiltonian(ii/countDiag))]);
        figure(1)
        plot(x,u); ylim([-.5 2.5]);
        title('solution of KdV equation, non-symplectic')
    end
    
    M = 1./(1-.25*dt*Lk);
    % First stage ARK4
    k0 = -3*Dk.*u2k_dealiasing(uk,params);
    l0 = Lk.*uk;
    % Second stage
    u1 = M.*(uk+.5*dt*k0+.25*dt*l0);
    k1 = -3*Dk.*u2k_dealiasing(u1,params);
    l1 = Lk.*u1;
    % Third stage
    u2 = M.*(uk+dt*(13861*k0/62500+6889*k1/62500+8611*l0/62500-1743*l1/31250));
    k2 = -3*Dk.*u2k_dealiasing(u2,params);
    l2 = Lk.*u2;
    % Fourth stage
    u3 = M.*(uk+dt*(-0.04884659515311858*k0-0.1777206523264010*k1+0.8465672474795196*k2...
    +0.1446368660269822*l0-0.2239319076133447*l1+0.4492950415863626*l2));
    k3 = -3*Dk.*u2k_dealiasing(u3,params);
    l3 = Lk.*u3;
    % Fifth stage
    u4 = M.*(uk+dt*(-0.1554168584249155*k0-0.3567050098221991*k1+1.058725879868443*k2...
    +0.3033959883786719*k3+0.09825878328356477*l0-0.5915442428196704*l1...
    +0.8101210538282996*l2+0.2831644057078060*l3));
    k4 = -3*Dk.*u2k_dealiasing(u4,params);
    l4 = Lk.*u4;
    % Sixth stage
    u5 = M.*(uk+dt*(0.2014243506726763*k0+0.008742057842904184*k1+0.1599399570716811*k2...
    +0.4038290605220775*k3+0.2260645738906608*k4+0.1579162951616714*l0...
    +0.1867589405240008*l2+0.6805652953093346*l3-0.2752405309950067*l4));
    k5 = -3*Dk.*u2k_dealiasing(u5,params);
    l5 = Lk.*u5;
    
    % Error control
%     r1 = dt*max(max(max(abs(ifft(0.003204494398459*(k0+l0) -0.002446251136679*(k2+l2)-0.021480075919587*(k3+l3)...
%         +0.043946868068572*(k4+l4) -0.023225035410765*(k5+l5))))));
%     if r1>tol,dt=.75*dt;continue,end
    
    % Successful step, proceed to evaluation
    t = t+dt;
    u = real(ifft(uk+dt*(0.1579162951616714*(k0+l0)+0.1867589405240008*(k2+l2)+...
    0.6805652953093346*(k3+l3)-0.2752405309950067*(k4+l4)+(k5+l5)/4)));
    uk = fft(u);

    % step size adjustment: EPS, PI.3.4 ; divide by 4 for a 4th order
    % method with 3rd order embedded
    % the adaptive time step is turned off
%     dt = ((.75*tol/r1)^(.3/4))*((r0/r1)^(.4/4))*dt;
%     r0=r1;
    if any(abs(u(:))>ulim),break,end
end


% plot results
figure(1);
plot(x,u,'LineWidth',2); hold on;
plot(x,2*sech(mod(x-4*t+Ls,2*Ls)-Ls).^2,'k--');
ylim([-.5 2.5]);
legend('numerics','theory')
title('solution of KdV equation, non-symplectic');

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