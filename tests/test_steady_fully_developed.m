% test fully developed
addpath(genpath('../source'));
[A, Fp, op] = discretize(100, 0, 4);
Dr2 = op.Dr2;
rm = op.rm;
f    = -1;
F   = f*ones(size(Dr2, 1), 1);
v   = Dr2\F;
va  = f/4*(rm.^2 - 1);
plot(rm, v, '-',rm, va,'--','linew',2);
set(gca,'fontsize', 16);
title('test steady fully developed flow')
legend('numerical','analytical')
xlabel('r');
ylabel('v');
shg



