function u2k = u2k_dealiasing(uk,p)
% Function takes Fourier coefficients of u and struct containing
% parameters (p) and evaluates the term u^2. 
% Returns Fourier coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% persistent dX DX
% if isempty(DX)
%     k = (pi/p.Ls)*[0:p.J/2-1 0 -p.J/2+1:-1]';
%     dX = 1i*k;
%     % For the dealiased jacobian:
%     k = (pi/p.Ls)*[0:.75*p.J-1 0 -.75*p.J+1:-1]';
%     DX = 1i*k;
%     clear k
% end

    
% For using a 3/2-rule dealiased jacobian:
    % physical space, 3/2 grid; factor of (3/2) scales fft
    U_hat = zeros([1.5*p.J 1]);
    U_hat(1:p.J/2) = (3/2)*uk(1:p.J/2);
    U_hat(p.J+2:1.5*p.J) = (3/2)*uk(p.J/2+2:p.J);
    % calculate u.gradq on 3/2 grid
    u = real(ifft(U_hat));
    u2 = u.^2;
    % fft, 3/2 grid; factor of (2/3) scales fft
    U2_hat = (2/3)*fft(u2);
    % reduce to normal grid
    u2k = zeros([p.J 1]);
    u2k(1:p.J/2) = U2_hat(1:p.J/2);
    u2k(p.J/2+2:p.J) = U2_hat(p.J+2:1.5*p.J);