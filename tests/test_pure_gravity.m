% test purely gravity mode

% test fully developed
addpath(genpath('../source'));
[A, Fp, op] = discretize(100, 0, 8);

% search the mode near the purely gravity mode.
[v, d] = eigs(A, 1, 1i*1);

fprintf('Pure gravity mode omega = %8.5f\n', imag(d));
plot(op.rm, real(v(1:end-1))./ max(abs(v(1:end-1))));