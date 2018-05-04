function [v, u] = solution_bessel(r, Lambda, Omega, F, m)
%  [v, u] = solution_bessel(r, Lambda, Omega, F, m)
% r : radial coordinates, vecter of [nr, 1];
% Lambda: non-dimensional parameter Lambda
% Omega: non-dimensional frequency.
% F: magnitude of the forcing at frequency Omega.
% m: the number of terms to be included in series expansion of bessel
% function.
%

alpha = 1i*Lambda/Omega;
M = 0:m;

% coefficient phi, [1, m]
phi_m = 0.5*(-1).^M./(factorial(M+1).*gamma(M+2)).*(-1i/2/sqrt(alpha)).^(2*M+1);

B = 1i*F/Omega/sum(phi_m.*(1 + (1 ./ (M + 2) - 1)*1/Omega^2));
x = -1i*r/sqrt(alpha);
u = B * besselj(1, x);

% note matlab broadcast is used here:
v = B * (r.^(2*M+2) - 1)*phi_m';
end

