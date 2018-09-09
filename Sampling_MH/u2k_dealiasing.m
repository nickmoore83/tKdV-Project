function u2k = u2k_dealiasing(uk)
% Function takes Fourier coefficients of u and struct containing
% parameters (p) and evaluates the term u^2. 
% Returns Fourier coefficients.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

J = length(uk);
   
% For using a 3/2-rule dealiased jacobian:
    % physical space, 3/2 grid; factor of (3/2) scales fft
    U_hat = zeros([1.5*J 1]);
    U_hat(1:J/2) = (3/2)*uk(1:J/2);
    U_hat(J+2:1.5*J) = (3/2)*uk(J/2+2:J);
    % calculate u.gradq on 3/2 grid
    u = real(ifft(U_hat));
    u2 = u.^2;
    % fft, 3/2 grid; factor of (2/3) scales fft
    U2_hat = (2/3)*fft(u2);
    % reduce to normal grid
    u2k = zeros([J 1]);
    u2k(1:J/2) = U2_hat(1:J/2);
    u2k(J/2+2:J) = U2_hat(J+2:1.5*J);