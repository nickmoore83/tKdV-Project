% Sampling from microcanonical-canonical ensemble with constant energy
clear;

Lambda = 20;  % number of spectral modes
E0 = 1;        % total energy
theta = -.5;  % inverse temperatur
D0 = .6;       % depth
Ls = pi;    % half domain size
Dk = 1i*(2*pi)/(2*Ls)*[0:Lambda/2 -Lambda/2+1:-1]'; % wavenumbers

MC = 1E7;     % starting # of samples
samp_H = zeros(MC,1);

tic;
for ii=1:MC
    uk = fft(randn(Lambda,1));
    E = .5*sum(abs(uk).^2)/Lambda^2;
    uk = sqrt(E0/E)*uk;
    
    samp_H(ii) = ( 1/6*D0^(-13/4)*E0^(1/2)*real(sum(u2k_dealiasing(uk).*conj(uk))) ...
                  -1/2*D0^(3/2)*sum(abs(Dk.*uk).^2) )/Lambda^2;
end
toc;

Hmax = max(exp(-theta*samp_H));
gibbs = exp(-theta*samp_H)/Hmax;
samp_eff=samp_H(gibbs>=rand(MC,1));
MC_eff=length(samp_eff);

figure(1)
[counts,centers] = hist(samp_eff,100);
plot(centers,counts/MC_eff/(centers(2)-centers(1))); hold on;
ylabel('PDF of Hamiltonian')
title(['Lambda = ',num2str(Lambda),', E0 = ',num2str(E0),', D0 = ', num2str(D0), ', N = ',num2str(MC_eff)]);