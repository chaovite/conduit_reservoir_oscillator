function [A, Fp, op] = discretize(nr, lambda, order)
% u = discretize(nr, lambda)
% discretize the non-dimensionalised gravity-inertia equations with domain
% r in [0, 1].
%
% nr: number grid points in r direction
% lambda: the dimensionless parameter that controls the viscous
% dissipation.
%

% non-dimensional radius.
R = 1; 
dr = R/nr;
truncate = 1; % strong boundary condition by truncating.
[rp,rm,Pp,Pm,Qp,Qm] = sbp_staggered_strong(order,nr, dr, truncate);
Rp = spdiags(rp,0,nr+1,nr+1);
Rm = spdiags(rm,0,nr,nr);

% Assemble second derivative in r direction.
Dr2  = inv(Rm)*inv(Pm)*Qm*Rp*inv(Pp)*Qp;
% first derivative in r direction 
Dr1 =  Pp\Qp; % 1D. [nr+1, nr]
er = ones([nr,1]);

% width averaging operators for v. v is on the staggered grid.
W1 = (2/R^2)*er'*Rm*Pm; % Cylindrical width averaged operator that operators on just one row of v, dimension [1, nr]
% output all the operators for discretization.
op.Dr2 = Dr2;
op.Dr1 = Dr1;
op.W1  = W1;
op.rm   = rm;
op.rp    = rp;
op.Pp   = Pp;
op.Pm  = Pm;

% construct the operators A and Fp. dq = A*q + Fp.
% q = [v, h].

A11 = lambda*Dr2;
A12 = -er;
A21 = W1;
A22 = sparse(1, 1);

Fp   = [-er; sparse(1,1)];
A     = [A11, A12; ...
           A21, A22];
end