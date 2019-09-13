function SolverKdV_SymplecticM4a_MC(MC, um,Nw, theta,gibd, fi)
% This script solves the rescaled KdV equation in a periodic domain
% using pseudo-spectral & symplectic M4a integrator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set simulation parameters
J = 32;      % Number of points in each direction, J devidable by 4
% MC = 1E4;  % reference MC size
dt = .5E-3;    % initial time step size
Nt = 1e5;    % Number of time steps
ulim = 1.5E5; % if any u > ulim, simulation stops

tol = 1e-12;  % iteration tolerance, ref 1e-10 1e-11 1e-20
w1=(2+2^(1/3)+2^(-1/3))/3;
w2=1-2*w1;     % integration fractional time steps

% model parameters from experiments
epsi0 = 0.016;  % amplitude-to-depth ratio
del0 = 0.22;    % depth-to-wavelength ratio
Dm = 1;  Dp = 0.24;  % depth ratio
if fi == 0     % fi == 0 for upstream flow
    D0 = Dm;
else
    D0 = Dp;
end
C2 = pi^2*del0/(6*Nw^2);
C3 = 3/2*pi^(1/2)*epsi0/del0;


% Put useful stuff into a struct
params = struct('MC',MC,'theta',theta, 'epsi0',epsi0,'del0',del0,'D0',D0,'Nw',Nw,...
                'C2',C2,'C3',C3, 'gibd',gibd,'J',J,'dt',dt,'um',um, 'Nt',Nt,'tol',tol);
x=(-pi:2*pi/J:pi-2*pi/J)';

% Set up implicit dispersion
k = repmat([0:J/2 -J/2+1:-1]',[1 MC]); % wavenumbers


% Initialize
[u0,uk0, enek,samp_H] = sampling_MH_phy(J,MC, theta,params);

t = 0;
uk0(1,:) = um*J;
uk = uk0;
ukm1 = uk0;  % at time n-1
ukm2 = uk0;  % at time n-2
% interpolatin function for the initial iteration input
interp = @(u0,um1,um2,w) (.5*w*(w+3)+1).*u0 - w*(w+2).*um1 +.5*w*(w+1).*um2;



% Diagnostics
countDiag = .01/dt;  % Compute diagnostics every countDiag steps
T = zeros(1,Nt/countDiag);
u = ifft(uk);
mass = zeros(1,Nt/countDiag);
energy = zeros(1,Nt/countDiag);
hamiltonian = zeros(1,Nt/countDiag);

meanu = zeros(J,Nt/countDiag);
covau = zeros(J,Nt/countDiag);
skewu = zeros(J,Nt/countDiag);
kurtu = zeros(J,Nt/countDiag);
skewdu = zeros(J,Nt/countDiag);
kurtdu = zeros(J,Nt/countDiag);



% Main loop 
% tic;
for ii=1:Nt
    if mod(ii,countDiag)==0
        if any(isnan(uk(:))), break, end
        T(ii/countDiag)=t;
        du=real(ifft(1i*k.*uk));
        
        meanu(:,ii/countDiag)=mean(uk,2);
        covau(:,ii/countDiag)=var(uk,[],2);
        skewu(:,ii/countDiag)=skewness(u,[],2);
        kurtu(:,ii/countDiag)=kurtosis(u,[],2);
        skewdu(:,ii/countDiag)=skewness(du,[],2);
        kurtdu(:,ii/countDiag)=kurtosis(du,[],2);
        
        mass(ii/countDiag) = mean( sum(u)/J );
        energy(ii/countDiag) = mean( 2*pi*.5*sum(abs(uk(2:end,:)).^2)/J^2 );
        temp = -(C3*D0^(-1)*1/6*real(sum(u2k_dealiasing_MC(uk,params).*conj(uk))) ...
                -C2*D0^(1)*1/2*sum(abs(k.*uk).^2)) *2*pi/J^2;
        hamiltonian(1:10,ii/countDiag) = temp(1:1);

        if mod(ii,1e3)==0
            display(['iteration i = ', num2str(ii),'; E = ',num2str(energy(ii/countDiag)),', H = ',num2str(mean(hamiltonian(:,ii/countDiag)))]);
        end
        % toc;
        
        % SAVE DATA TO OUTPUT FILE
        save('output.mat', 'T','meanu','covau','skewu','kurtu', 'skewdu','kurtdu', 'energy','hamiltonian')
    end
    
    
    % First stage
    zk0 = .5*(uk+interp(uk,ukm1,ukm2,w1));  % initial iteration
    [zk,ik] = iterativeKdV_M4a_MC(zk0, w1*dt,uk,params, tol);
    if ik==1
        disp('Error! quit...');
        break;
    end
    yk = 2*zk-uk;
    % Second stage
    zk0 = .5*(yk+interp(uk,ukm1,ukm2,w1+w2));
    [zk,ik] = iterativeKdV_M4a_MC(zk0, w2*dt,yk,params, tol);
    if ik==1
        disp('Error! quit...');
        break;
    end
    yk = 2*zk-yk;
    % Third stage
    zk0=.5*(yk+interp(uk,ukm1,ukm2,1));
    [zk,ik] = iterativeKdV_M4a_MC(zk0, w1*dt,yk,params, tol);
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
    uk = fft(u); uk(J/2+1,:)=0;

    if any(abs(u(:))>ulim),break,end
end

