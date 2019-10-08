% Sampling from microcanonical-canonical ensemble with constant energy
% using a Metropolis-Hastings MC approach
function [samp_H0,samp_H1,samp_u,enek] = sampling_MH_scaled(D0,D1,Nw, theta, MC)
% Sampling from microcanonical-canonical ensemble with constant energy
% using a Metropolis-Hastings MC approach
J=32; %Nw=10; theta=-30.; D0=.24; 
%MC = 1E4;
% model parameters from experiments
epsi0 = 0.017;  % amplitude-to-depth ratio
del0 = 0.22;    % depth-to-wavelength ratio
C2 = (2/3)*pi^2*del0/Nw^2;
C3 = (3/2)*pi^(1/2)*epsi0/del0;


kk = [0:J/2 -J/2+1:-1]';
Dk = 1i*kk; % wavenumbers
Enorm = 1;      % unit energy in modes

count = 10;
cstart = 10000;
N_iter = cstart+MC*count;    % sampling steps along the Markov chain
samp_H0 = zeros(1,MC);        % sampled Hamiltonian
samp_H1 = zeros(1,MC); 
samp_u = zeros(J,MC);        % sampled states
enek = zeros(J,1);  % energy spectrum for each mode


% initial setup of the sample
uk = fft(randn(J,1));  uk(1)=0; uk(J/2+1)=0;
E = .5*sum(abs(uk).^2) *2*pi/J^2;
uk = sqrt(Enorm/E)*uk; 

% previous distribution
H_pre = ( 1/6*D0^(-3/2)*C3*real(sum(u2k_dealiasing(uk).*conj(uk))) ...
          -1/2*D0^(1/2)*C2*sum(abs(Dk.*uk).^2) ) *2*pi/J^2;
count_acpt = zeros(N_iter,1);
% run Markov chain p(x|y)~exp(-.5*|x-y|^2)
ap = 1.;
% tic;
for ii = 1:N_iter

    % update componentwisely
    for jj = 1:J
        u_tilde = ifft(uk);

        pert = ap*sqrt(Enorm/J)*randn(1);
        u_tilde(jj) = u_tilde(jj)+pert;
        uk_tilde = fft(u_tilde); uk_tilde(1)=0; uk_tilde(J/2+1)=0;
        E = .5*sum(abs(uk_tilde).^2) *2*pi/J^2;
        uk_tilde = sqrt(Enorm/E)*uk_tilde;

        
        % new distribution
        H_new = ( 1/6*D0^(-3/2)*C3*real(sum(u2k_dealiasing(uk_tilde).*conj(uk_tilde))) ...
                  -1/2*D0^(1/2)*C2*sum(abs(Dk.*uk_tilde).^2) ) *2*pi/J^2;
        alpha = exp(-theta*(H_new-H_pre));
        r = rand(1);
        if r <= alpha
            uk = uk_tilde;
            H_pre = H_new;
            count_acpt(ii) = count_acpt(ii)+1;
        end
    end
    
    % record sample
    if ii > cstart && mod(ii-cstart,count)==0
        samp_H0((ii-cstart)/count) = H_pre;
        samp_u(:,(ii-cstart)/count) = ifft(uk);
        enek = enek+ .5*abs(uk).^2 *2*pi/J^2;
        
        samp_H1((ii-cstart)/count) = ( 1/6*D1^(-3/2)*C3*real(sum(u2k_dealiasing(uk_tilde).*conj(uk_tilde))) ...
                  -1/2*D1^(1/2)*C2*sum(abs(Dk.*uk_tilde).^2) ) *2*pi/J^2;
    end
end
% toc;
enek = enek/MC;