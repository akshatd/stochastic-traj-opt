function [Aext, Bext, Q, R, P] = extendState(system, Qe, Re, Pe)

A = system.A;
B = system.B;
C = system.C;

nx = size(A,1);  % Number of states
nu = size(B,2);  % Number of inputs
ny = size(C,1);  % Number of measurements
nr = ny; % Number of references (equivalent to measurments)

% Extend states are [x_k, u_k-1, r_k]
Aext = [A,            B,       zeros(nx, nr);
	zeros(nu,nx), eye(nu), zeros(nu, nr);
	zeros(nr,nx), zeros(nr,nu), eye(nr)];

Bext = [B;
	eye(nu, nu);
	zeros(nr,nu)];

E = [C, zeros(ny,nu), -eye(nr)];

Q = E'*Qe*E;
R = Re;
P = blkdiag(Pe, zeros(size(Aext,1) - nx));

end
