% study parametric study of gravity oscillation mode
source = '../source/gravity_inertia_oscillator/source/';
addpath(genpath(source));
%% 
tic
n  = 200;
chis = logspace(0, -4, n); % the parameter \chi in the paper part 1.

% Ts in this script is non-dimensional T*= T/T0 in the paper part 1.
 
nr      = 100;
order = 8;
Ts = zeros(1, n);
Qs          = zeros(1, n);
for i = 1: n
    [A, Fp, op] = discretize(nr, chis(i), order);
    omega = eigs(A, 1, 1i);
    T = 1/imag(omega);
    Q  = abs(imag(omega)/real(omega))/2;
    Ts(i) = T;
    Qs(i) = Q;
end
toc
%% plot the nondimensional period and quality factor v.s. dimensionless parameter Chi.
figure(1);
ind = Qs>0.5;
Chi_critical = chis(Qs<=0.5);
Chi_critical = Chi_critical(end);

yyaxis right
semilogx(chis(ind), Ts(ind),'linew',2);
ylabel('$T^*=T/T_0$', 'Interpreter','latex');
xlabel('$\chi$', 'Interpreter','latex');
set(gca,'fontsize',20);

yyaxis left
semilogx(chis(ind), Qs(ind),'linew',2);
hold on;
semilogx(Chi_critical, 0, 'kd','markersize',15,'markerfacecolor','k');
ylabel('Quality factor $Q$', 'Interpreter','latex');
grid minor;
hold off