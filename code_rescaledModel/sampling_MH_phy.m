function [samp_u,sampk_u, enek,samp_H] = sampling_MH_phy(J, MC, theta,p)
% Sampling from microcanonical-canonical ensemble with constant energy
% using a Metropolis-Hastings MC approach
C2=p.C2; C3=p.C3;
D0=p.gibd;

kk = [0:J/2 -J/2+1:-1]';
Dk = 1i*kk; % wavenumbers
Enorm = 1;      % unit energy in modes

count = 10;
cstart = 10000;
N_iter = cstart+MC*count;    % sampling steps along the Markov chain
samp_H = zeros(1,MC);        % sampled Hamiltonian
samp_u = zeros(J,MC);        % sampled states
sampk_u = zeros(J,MC);
enek = zeros(J,1);  % energy spectrum for each mode


% initial setup of the sample
uk = fft(randn(J,1));  uk(1)=0; uk(J/2+1)=0;
E = .5*sum(abs(uk).^2) *2*pi/J^2;
uk = sqrt(Enorm/E)*uk; 

% previous distribution
H_pre = ( 1/6*D0^(-1)*C3*real(sum(u2k_dealiasing(uk,p).*conj(uk))) ...
          -1/2*D0^(1)*C2*sum(abs(Dk.*uk).^2) ) *2*pi/J^2;
count_acpt = zeros(N_iter,1);
% run Markov chain p(x|y)~exp(-.5*|x-y|^2)
ap = 1.;
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
        H_new = ( 1/6*D0^(-1)*C3*real(sum(u2k_dealiasing(uk_tilde,p).*conj(uk_tilde))) ...
                  -1/2*D0^(1)*C2*sum(abs(Dk.*uk_tilde).^2) ) *2*pi/J^2;
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
        samp_H((ii-cstart)/count) = H_pre;
        samp_u(:,(ii-cstart)/count) = ifft(uk);
        sampk_u(:,(ii-cstart)/count) = uk;
        enek = enek+ .5*abs(uk).^2 *2*pi/J^2;
    end
end
enek = enek/MC;