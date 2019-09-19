function ham = hamiltonian(uk,params)
% Basic stuff.
J = size(uk,1);
MC = size(uk,2);
kvec = [0:J/2 -J/2+1:-1]';
kmat = repmat(kvec,[1 MC]);
dx = 2*pi/J;
% Compute the Hamiltonian piece by piece.
ham2 = 0.5 * sum(abs(kmat.*uk).^2) * dx;
ham3 = 1/6 * real(sum(u2k_dealiasing_MC(uk,params).*conj(uk))) * dx;
ham = C2*D0^(1/2)*ham2 - C3*D0^(-3/2)*ham3;


% SolverKdV_SymplecxticM4a_MC.m lines 102-103
%k = repmat([0:J/2 -J/2+1:-1]',[1 MC]); % wavenumbers
%temp = -(C3*D0^(-1)*1/6*real(sum(u2k_dealiasing_MC(uk,params).*conj(uk))) ...
 %       -C2*D0^(1)*1/2*sum(abs(k.*uk).^2)) *2*pi/J^2;

 % sampling_MH_phy.m lines 26-27
%kk = [0:J/2 -J/2+1:-1]';
%Dk = 1i*kk; % wavenumbers
%H_pre = -( 1/6*D0^(-1)*C3*real(sum(u2k_dealiasing(uk,p).*conj(uk))) ...
%      -1/2*D0^(1)*C2*sum(abs(Dk.*uk).^2) ) *2*pi/J^2;

% sampling_MH_phy.m lines 45-46
%H_new = -( 1/6*D0^(-1)*C3*real(sum(u2k_dealiasing(uk_tilde,p).*conj(uk_tilde))) ...
%  -1/2*D0^(1)*C2*sum(abs(Dk.*uk_tilde).^2) ) *2*pi/J^2;
