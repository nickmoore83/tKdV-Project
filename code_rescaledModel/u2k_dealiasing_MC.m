function u2k = u2k_dealiasing_MC(uk)
% Function takes Fourier coefficients of u and struct containing
% parameters (p) and evaluates the term u^2. 
% Returns Fourier coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Modified by Nick Moore on Sep 19, 2016.
%% Key: if J and MC are defined by the size of u, then this code works for 
%% both u2k_dealiasing_MC and for the regular u2k_dealiasing.

% persistent dX DX
% if isempty(DX)
%     k = (pi/p.Ls)*[0:JJ/2-1 0 -JJ/2+1:-1]';
%     dX = 1i*k;
%     % For the dealiased jacobian:
%     k = (pi/p.Ls)*[0:.75*JJ-1 0 -.75*JJ+1:-1]';
%     DX = 1i*k;
%     clear k
% end

    
% For using a 3/2-rule dealiased jacobian:
    % physical space, 3/2 grid; factor of (3/2) scales fft
    JJ = size(uk,1);
    MC = size(uk,2);
    U_hat = zeros([1.5*JJ MC]);    
    U_hat(1:JJ/2,:) = (3/2)*uk(1:JJ/2,:);
    U_hat(JJ+2:1.5*JJ,:) = (3/2)*uk(JJ/2+2:JJ,:);
    % calculate u.gradq on 3/2 grid
    u = real(ifft(U_hat));
    u2 = u.^2;
    % fft, 3/2 grid; factor of (2/3) scales fft
    U2_hat = (2/3)*fft(u2);
    % reduce to normal grid
    u2k = zeros([JJ MC]);
    u2k(1:JJ/2,:) = U2_hat(1:JJ/2,:);
    u2k(JJ/2+2:JJ,:) = U2_hat(JJ+2:1.5*JJ,:);