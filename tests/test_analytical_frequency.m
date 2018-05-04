% compare with analytical solution in frequency domain.
tic
% add source;
addpath(genpath('../source'));

% parameters
Lambda = 1e-4;
Omega  = 1; % the frequency of gravity-ineria mode.
fhat  = 1; % source amplitude.

%% numerical solution
nr           = 200;
order      = 4;
tic
[A, Fp, op] = discretize(nr, Lambda, order);
B      = (-1i*Omega*speye(size(A,1)) - A);
toc
qhat = B\(Fp*fhat);
vhat = qhat(1:end-1);
uhat = op.Dr1*vhat;
uhat(1) = 0;

% figure(1)
% plot(op.rm, abs(vhat));
% 
% figure(2)
% plot(op.rp, abs(uhat));
% analytical solutions
r = linspace(0, 1, 201)';
[v, u] = solution_bessel(r, Lambda, Omega, fhat, 95);
%%
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
plot(r, abs(un),'-', op.rp, abs(uhat_n),'--','linew',1.5);
%
figure(2)
plot(r, abs(vn),'-', op.rm, abs(vhat_n),'--','linew',1.5);
%% plot the phase difference:
% figure(3)
% plot(r, abs(v./uhat),'-','linew',1.5);
% ylim([-10, 10])
