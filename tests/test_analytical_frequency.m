% compare with analytical solution in frequency domain.
tic
% add source;
addpath(genpath('../source'));

% parameters
Lambda = 1e-2;
Omega  = 1; % the frequency of gravity-ineria mode.
fhat  = 1; % source amplitude.

% numerical solution
nr           = 500;
order      = 4;
tic
[A, Fp, op] = discretize(nr, Lambda, order);
B      = (-1i*Omega*speye(size(A,1)) - A);
toc
qhat = B\(Fp*fhat);
vhat = qhat(1:end-1);
uhat = op.Dr1*vhat;
uhat(1) = 0;

[~, u] = solution_bessel(op.rp, Lambda, Omega, fhat, 90);
[v, ~] = solution_bessel(op.rm, Lambda, Omega, fhat, 90);
%
vhat_n = vhat./max(abs(vhat));
uhat_n = uhat./max(abs(uhat));
vn        = v/max(abs(v));
un        = u/max(abs(u));

% compare v
% plot(r, real(vn), op.rm, real(vhat_n), r, imag(vn), op.rm, imag(vhat_n));
% legend()

% compare u
% plot(r, real(un), op.rp, real(uhat_n), r, imag(un), op.rp, imag(uhat_n));
figure(1)
plot(op.rp, abs(un),'-', op.rp, abs(uhat_n),'--','linew',1.5);
%
figure(2)
plot(op.rm, abs(vn),'-', op.rm, abs(vhat_n),'--','linew',1.5);
%% plot the phase difference:
figure(1);
plot(op.rp, imag(u),'-', op.rp, imag(uhat),'--','linew',1.5);
figure(2);
plot(op.rp, real(u),'-', op.rp,real(uhat),'--','linew',1.5);
%%
figure(3);
plot(op.rm, imag(v),'-', op.rm, imag(vhat),'--','linew',1.5);
figure(4);
plot(op.rm, real(v),'-', op.rm,real(vhat),'--','linew',1.5);
