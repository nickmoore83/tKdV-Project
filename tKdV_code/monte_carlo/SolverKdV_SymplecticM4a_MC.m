function SolverKdV_SymplecticM4a_MC(J, Ls,D0,E0, um, MC, theta,gibd)
% This script solves KdV equation in a periodic domain
% using pseudo-spectral & symplectic M4a integrator
% referecen parameters;
% J = 32; Ls = 3; D0 = 1, 0.24; E0 = 50, 100, 150; um = 0;
% MC = 1e4; theta = -.5; gibd = 1, 0.24;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set simulation parameters
dt = .5E-3 *J/32;    % initial time step size
Nt = 10e5 *64/J;    % Number of time steps
ulim = 1.5E5; % if any u > ulim, simulation stops

tol = 1e-12;  % iteration tolerance, ref 1e-10 1e-11 1e-20
w1=(2+2^(1/3)+2^(-1/3))/3;
w2=1-2*w1;     % integration fractional time steps

alpha1=-3/2; alpha2=1/2;


% Put useful stuff into a struct
params = struct('MC',MC,'theta',theta, 'Ls',Ls,'D0',D0,'E0',E0,'alpha1',alpha1,'alpha2',alpha2,...
                'gibd',gibd,'J',J,'dt',dt,'um',um, 'Nt',Nt,'tol',tol);
x=(-pi*Ls:2*pi*Ls/J:pi*Ls-2*pi*Ls/J)';

% Set up implicit dispersion
k = repmat([0:J/2 -J/2+1:-1]',[1 MC]); % wavenumbers


% Initialize
[u0,uk0, enek,samp_H] = sampling_MH_phy(J,Ls,gibd,E0, um, MC, theta,params);

t = 0;
% u0 = randn(J,MC); 
% uk0 = fft(u0);  uk0(1,:)=0; uk0(J/2+1,:)=0;
% ene0=ones(J,1)*sqrt(2*pi*.5*sum(abs(uk0).^2)) /J;
% uk0 = uk0./ene0;
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
hamiltonian = zeros(10,Nt/countDiag);

meanu = zeros(J,Nt/countDiag);
covau = zeros(J,Nt/countDiag);
skewu = zeros(J,Nt/countDiag);
kurtu = zeros(J,Nt/countDiag);

% save data in binary file
% fid=fopen(['/scratch/dq271/tKdV_simu/dist/tKdV_MC',num2str(MC/100),'th',num2str(abs(theta*10)),'gd',num2str(gibd*100),'_dt',num2str(dt*1e5),'J',num2str(J),'Ls',num2str(Ls*10),'E',num2str(E0*10),'D',num2str(D0*100),'Um',num2str(um*100),'.bin'],'w');
% fwrite(fid,u,'double');
% fclose(fid);


% Main loop 
% tic;
for ii=1:Nt
    if mod(ii,countDiag)==0
        if any(isnan(uk(:))), break, end
        T(ii/countDiag)=t;
        
        meanu(:,ii/countDiag)=mean(uk,2);
        covau(:,ii/countDiag)=var(uk,[],2);
        skewu(:,ii/countDiag)=skewness(u,[],2);
        kurtu(:,ii/countDiag)=kurtosis(u,[],2);
        
        mass(ii/countDiag) = mean( sum(u)/J );
        energy(ii/countDiag) = mean( 2*pi*.5*sum(abs(uk(2:end,:)).^2)/J^2 );
        temp = (E0^(1/2)*(2*Ls)^(-3/2)*D0^(alpha1)*1/6*real(sum(u2k_dealiasing_MC(uk,params).*conj(uk))) ...
                                     -(2*Ls)^(-3)*D0^(alpha2)*1/2*sum(abs(k.*uk).^2)) *2*pi/J^2;
        hamiltonian(1:10,ii/countDiag) = temp(1:10);

%         fid=fopen(['/scratch/dq271/tKdV_simu/dist/tKdV_MC',num2str(MC/100),'th',num2str(abs(theta*10)),'gd',num2str(gibd*100),'_dt',num2str(dt*1e5),'J',num2str(J),'Ls',num2str(Ls*10),'E',num2str(E0*10),'D',num2str(D0*100),'Um',num2str(um*100),'.bin'],'a');
%         fwrite(fid,u,'double');
%         fclose(fid);
        if mod(ii,10e3 *64/J)==0
%             save(['/scratch/dq271/tKdV_simu/dist/tKdV_MC',num2str(MC/100),'th',num2str(abs(theta*10)),'gd',num2str(gibd*100),'_dt',num2str(dt*1e5),'J',num2str(J),'Ls',num2str(Ls*10),'E',num2str(E0*10),'D',num2str(D0*100),'Um',num2str(um*100),'.mat'],...
%                  'ii','countDiag','dt','params','x','T','k', 'mass','energy','hamiltonian','uk','u', ...
%                  'meanu','covau','skewu','kurtu', 'u0','enek','samp_H');
            display(['iteration i = ', num2str(ii),'; E = ',num2str(energy(ii/countDiag)),', H = ',num2str(mean(hamiltonian(:,ii/countDiag)))]);
        end
%             toc;
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



% plot results
% figure(1)
% plot(T,mass); hold on;
% title('time-series of the total mass');
% figure(2)
% plot(T,energy); hold on;
% title('time-series of the total energy');
% figure(3)
% plot(T,hamiltonian); hold on;
% title('time-series of the Hamiltonian');
% 
% figure
% contour(x,T,ut',50); colorbar;
% xlabel('x'); ylabel('t');
% title('time-series of the KdV solution')
% 
% figure
% plot(0:Dt:Dt*Tlag,corp_norm(:,1),'-'); hold on
% plot(0:Dt:Dt*Tlag,corp_norm(:,2),'-');
% plot(0:Dt:Dt*Tlag,corp_norm(:,3),'-');
% plot(0:Dt:Dt*Tlag,corp_norm(:,4),'-');
% figure
% plot(0:Dt:Dt*Tlag,cors_norm(:,1),'-'); hold on
% plot(0:Dt:Dt*Tlag,cors_norm(:,2),'-');
% plot(0:Dt:Dt*Tlag,cors_norm(:,3),'-');
% plot(0:Dt:Dt*Tlag,cors_norm(:,4),'-');
% figure
% plot(linspace(-pi,pi,J),varp,'.-')
% figure
% plot(k(2:J/2+1),vars,'o-')
