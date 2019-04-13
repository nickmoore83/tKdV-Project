% This script solves KdV equation in a periodic domain
% u_t=-uu_x-u_xxx, -L<=x<=L (with model coefficients)
% using pseudo-spectral & symplectic M4a integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model parameters
D0 = .24;        % 1, 0.24, water depth
E0 = 100;       % 50, 200, model total energy
Ls = 3;        % half domain size
um = 0;        % initial mean state

% Set simulation parameters
J = 32;        % Number of points in each direction, J devidable by 4
dt = 1E-3;     % initial time step size
Nt = 50e3;     % Number of time steps
ulim = 1.5E5;  % if any u > ulim, simulation stops


% parameters for symplectic integrator
tol = 1e-12;   % iteration tolerance
w1=(2+2^(1/3)+2^(-1/3))/3;
w2=1-2*w1;     % integration fractional time steps
alpha1=-3/2; alpha2=1/2;

% Put useful stuff into a struct
params = struct('Ls',Ls,'D0',D0,'E0',E0,'alpha1',alpha1,'alpha2',alpha2, 'J',J,'dt',dt,'um',um, 'Nt',Nt,'tol',tol);
x=(-pi*Ls:2*pi*Ls/J:pi*Ls-2*pi*Ls/J)';

% Set up implicit dispersion
k = [0:J/2 -J/2+1:-1]'; % wavenumbers

% Initialize
theta=-.5;
[u0,uk0, enek,samp_H] = sampling_MH_phy(J,Ls,D0,E0, um, 1e3, theta,params);
t = 0; 
uk0 = uk0(:,end);
uk0(1)=0; uk0(J/2+1)=0;
uk0 = uk0/sqrt(2*pi*.5*sum(abs(uk0).^2)) *J;
uk0(1) = um*J;
uk = uk0;
ukm1 = uk0;  % at time n-1
ukm2 = uk0;  % at time n-2
% interpolatin function for the initial iteration input
interp = @(u0,um1,um2,w) (.5*w*(w+3)+1).*u0 - w*(w+2).*um1 +.5*w*(w+1).*um2;


% Diagnostics
countDiag = .01/dt;  % Compute diagnostics every countDiag steps
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
        T(ii/countDiag)=t;
        ut(:,ii/countDiag)=u;
        
        mass(ii/countDiag)=sum(u)/J;
        energy(ii/countDiag) = 2*pi*.5*sum(abs(uk(2:end)).^2)/J^2;
        hamiltonian(ii/countDiag) = (E0^(1/2)*(2*Ls)^(-3/2)*D0^(alpha1)*1/6*real(sum(u2k_dealiasing(uk,params).*conj(uk))) ...
                                     -(2*Ls)^(-3)*D0^(alpha2)*1/2*sum(abs(k.*uk).^2)) *2*pi/J^2;


%         toc
        if ~mod(ii,100)
        display(['iteration i = ', num2str(ii),'; E = ',num2str(energy(ii/countDiag)),...
                 ', H = ',num2str(hamiltonian(ii/countDiag))]);
        end
    end
    
    
    % First stage
    zk0 = .5*(uk+interp(uk,ukm1,ukm2,w1));  % initial iteration
    [zk,ik] = iterativeKdV_M4a(zk0, w1*dt,uk,params, tol);
    if ik==1
        disp('Error! quit...');
        break;
    end
    yk = 2*zk-uk;
    % Second stage
    zk0 = .5*(yk+interp(uk,ukm1,ukm2,w1+w2));
    [zk,ik] = iterativeKdV_M4a(zk0, w2*dt,yk,params, tol);
    if ik==1
        disp('Error! quit...');
        break;
    end
    yk = 2*zk-yk;
    % Third stage
    zk0=.5*(yk+interp(uk,ukm1,ukm2,1));
    [zk,ik] = iterativeKdV_M4a(zk0, w1*dt,yk,params, tol);
    if ik==1
        disp('Error! quit...');
        break;
    end
    yk = 2*zk-yk;

 
    ukm2 = ukm1;
    ukm1 = uk;
    % Successful step, proceed to evaluation
    uk = yk;
    t = t+dt;
    u = real(ifft(uk));
    uk = fft(u); uk(J/2+1)=0;

    if any(abs(u(:))>ulim),break,end
end



% plot results
figure
subplot(3,1,1)
plot(T,mass); hold on;
title('time-series of the total mass');
subplot(3,1,2)
plot(T,energy); hold on;
title('time-series of the total energy');
subplot(3,1,3)
plot(T,hamiltonian); hold on;
title('time-series of the Hamiltonian');

figure
plot(T,ut(end/2,:));
xlabel('time'); ylabel('u');
title('time series of u evaluated at x=0')

figure
contour(x,T,ut',50); colorbar;
xlabel('x'); ylabel('t');
title('time-series of the KdV solution')