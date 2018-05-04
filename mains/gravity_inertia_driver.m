% gravity-inertia oscillation driver
clear
addpath(genpath('../source'));
[A, Fp, op] = discretize(100, 1, 8);