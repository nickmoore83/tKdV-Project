function SolverKdV_SymplecticM4a_MC(C2,C3,Drat,Lambda,...
    MC,theta,gibd,fi,dt,nout,tfin)  

% This script solves the rescaled KdV equation in a periodic domain
% using pseudo-spectral & symplectic M4a integrator

% Input
% C2 and C3: Constants in Hamiltonian.
% Drat: Depth ratio; reference Drat = 0.24.
% Lambda: Number of modes; reference Lambda = 16.
% MC: The number of trajectories; reference MC = 1E4.
% theta: The inverse temperature; range 10-20.
% gibd: The Gibbs distribution to sample the initial state.
% fi: Flow index; fi = 0 for incoming or 1 for outgoing.
% dt: The numerical time step; reference dt = .5E-3.
% nout: The number of output times; reference nout = 100.
% tfin: The stopping time of the simulation; reference tfin = 10.

% Set simulation parameters
JJ = 2*Lambda;          % Number of grid points.
Nt = round(tfin/dt);    % Number of time steps.
MM = round(Nt/nout);    % Save output every MM steps.
um = 0.;                % zero momentum
% Other simulation parameters
ulim = 1.5E5; % if any u > ulim, simulation stops
tol = 1e-12;  % iteration tolerance, ref 1e-10 1e-11 1e-20
w1=(2+2^(1/3)+2^(-1/3))/3;
w2=1-2*w1;     % integration fractional time steps
kvec = repmat([0:JJ/2 -JJ/2+1:-1]',[1 MC]); % wavenumbers
% fi == 0 for upstream flow
if fi == 0     
    D0 = 1;
else
    D0 = Drat;
end
% Put useful stuff into a struct
params = struct('C2',C2,'C3',C3,'D0',D0,'J',JJ,'MC',MC,'gibd',gibd);
% Initialize
[u0,uk0, enek,samp_H] = sampling_MH_phy(JJ,MC, theta,params);
% Initialize
tt = 0.;
uk0(1,:) = um*JJ;
uk = uk0;
ukm1 = uk0;  % at time n-1
ukm2 = uk0;  % at time n-2
% interpolatin function for the initial iteration input
interp = @(u0,um1,um2,w) ...
    (.5*w*(w+3)+1).*u0 - w*(w+2).*um1 +.5*w*(w+1).*um2;
% Main variables
tsave = zeros(1,Nt/MM);
uu = ifft(uk);
uarray = zeros(JJ,MC,Nt/MM);
duarray = zeros(JJ,MC,Nt/MM);

% Main loop 
for ii=1:Nt
    if mod(ii,MM)==0
        if any(isnan(uk(:))), break, end
        iout = ii/MM;
        du=real(ifft(1i*kvec.*uk));
        % Save the microstate, u and du.
        tsave(iout)=tt;
        uarray(:,:,iout) = uu;
        duarray(:,:,iout) = du;
        % Check energy and hamiltonian.        
        energy = mean( 2*pi*.5*sum(abs(uk(2:end,:)).^2)/JJ^2 );
        ham = mean( hamiltonian(uk,C2,C3,D0) );
        % Print progress to screen.
        display(['step ',num2str(ii),...
            ', completed ',num2str(100*ii/Nt),'%']);
        display(['E = ',num2str(energy),', H = ',num2str(ham),newline]);
    end
    
    % First stage
    zk0 = .5*(uk+interp(uk,ukm1,ukm2,w1));  % initial iteration
    [zk,ik] = iterativeKdV_M4a_MC(zk0, w1*dt,uk,params, tol);
    if ik==1
        disp('Error! quit...'); break;
    end
    yk = 2*zk-uk;
    % Second stage
    zk0 = .5*(yk+interp(uk,ukm1,ukm2,w1+w2));
    [zk,ik] = iterativeKdV_M4a_MC(zk0, w2*dt,yk,params, tol);
    if ik==1
        disp('Error! quit...'); break;
    end
    yk = 2*zk-yk;
    % Third stage
    zk0=.5*(yk+interp(uk,ukm1,ukm2,1));
    [zk,ik] = iterativeKdV_M4a_MC(zk0, w1*dt,yk,params, tol);
    if ik==1
        disp('Error! quit...'); break;
    end
    yk = 2*zk-yk;
    ukm2 = ukm1;
    ukm1 = uk;
    % Successful step, proceed to evaluation
    uk = yk;
    tt = tt+dt;
    uu = real(ifft(uk));
    uk = fft(uu); uk(JJ/2+1,:)=0;
    if any(abs(uu(:))>ulim),break,end
end

% Save output
save('output.mat', 'C2','C3','Drat','Lambda',...
    'MC','theta','gibd','fi','dt','nout','tfin',...
    'uarray', 'duarray', 'tsave');