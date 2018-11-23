% Sampling from microcanonical-canonical ensemble with constant energy
% using a Metropolis-Hastings MC approach
clear;

Lambda = 20;  % number of spectral modes
E0 = 1;        % total energy
theta = -0.5;  % inverse temperatur
D0 = .6;       % depth
Ls = pi;    % half domain size
kk = (2*pi)/(2*Ls)*[0:Lambda/2 -Lambda/2+1:-1]';
Dk = 1i*kk; % wavenumbers

MC = 5E4;     % # of samples to pick
count = 10;
cstart = 10000;
N_iter = cstart+MC*count;    % sampling steps along the Markov chain
samp_H = zeros(MC,1);        % sampled Hamiltonian
enek = zeros(Lambda,1);  % energy spectrum for each mode


% initial setup of the sample
uk = fft(randn(Lambda,1));  %uk(Lambda/2+1)=0;
E = .5*sum(abs(uk).^2)/Lambda^2;
uk = sqrt(E0/E)*uk; 

% previous distribution
H_pre = ( 1/6*D0^(-13/4)*E0^(1/2)*real(sum(u2k_dealiasing(uk).*conj(uk))) ...
          -1/2*D0^(3/2)*sum(abs(Dk.*uk).^2) )/Lambda^2;
count_acpt = zeros(N_iter,1);
% run Markov chain p(x|y)~exp(-.5*|x-y|^2)
tic;
for ii = 1:N_iter

    % update componentwisely
    for jj = 1:Lambda/2+1
        uk_tilde = uk;
              
        if jj==1 || jj==Lambda/2+1
        pert = sqrt(.1*2*Lambda*E0)*randn(1);
        uk_tilde(jj) = uk(jj)+pert;
        E = .5*sum(abs(uk_tilde).^2)/Lambda^2;
        uk_tilde = sqrt(E0/E)*uk_tilde;
        else
        pert = sqrt(.1*2*Lambda*E0)*(randn(1)+1i*randn(1));
        uk_tilde(jj) = uk(jj)+pert;
        uk_tilde(end-jj+2) = conj(uk_tilde(jj));
        E = .5*sum(abs(uk_tilde).^2)/Lambda^2;
        uk_tilde = sqrt(E0/E)*uk_tilde;
        end
        
        % new distribution
        H_new = ( 1/6*D0^(-13/4)*E0^(1/2)*real(sum(u2k_dealiasing(uk_tilde).*conj(uk_tilde))) ...
                  -1/2*D0^(3/2)*sum(abs(Dk.*uk_tilde).^2) )/Lambda^2;
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
        enek = enek+abs(uk).^2/Lambda^2;
    end
end
toc;
enek = enek/MC;

figure(1)
[counts,centers] = hist(samp_H,100);
plot(centers,counts/MC/(centers(2)-centers(1))); hold on;
ylabel('PDF of Hamiltonian')
title(['Lambda = ',num2str(Lambda),', E0 = ',num2str(E0),', D0 = ', num2str(D0), ', N = ',num2str(MC)]);
figure(2)
plot(kk(1:Lambda/2+1),enek(1:Lambda/2+1),'-o'); hold on;
title('energy spectra <|u_k|^2>')