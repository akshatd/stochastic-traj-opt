% solve LQR problem and return optimal control gain
function Kopt = solveLQR(N, system, Q, R, P)
% Compute the matrices S, M, Qbar, Rbar, and K0N
% for the unconstrained LQ-MPC problem
%
% Inputs:
%   N: Prediction horizon
%   system: state space system
%   Q: State cost matrix
%   R: Input cost matrix
%   P: Terminal state cost matrix

A = system.A;
B = system.B;

nx = size(A, 1);
nu = size(B, 2);

% Initialize matrices
S = zeros(N * nx, N * nu);
M = zeros(N * nx, nx);
Qbar = zeros(N * nx, N * nx);
Rbar = zeros(N * nu, N * nu);

% Compute the first column of S
for i = 1:N
  rowStart = (i - 1) * nx + 1;
  rowEnd = i * nx;
  S(rowStart:rowEnd, 1:nu) = A^(i - 1) * B;
end

% Pad the first column and set it to other columns of S
for i = 2:N
  colStart = (i - 1) * nu + 1;
  colEnd = i * nu;
  zeroRows = (i - 1) * nx;
  zeroCols = nu;
  S(:, colStart:colEnd) = [zeros(zeroRows, zeroCols); S(1:end - zeroRows, 1:nu)];
end

% Compute first row of M
M(1:nx, :) = A;
% Compute the rest of M
for i = 2:N
  rowStart = (i - 1) * nx + 1;
  rowEnd = i * nx;
  % just multiply the previous rows by A to get higher powers
  M(rowStart:rowEnd, :) = A * M(rowStart - nx:rowEnd - nx, :);
end

% Compute Qbar except for the last row
for i = 1:N
  % Q is square so we can reuse indices
  rowStart = (i - 1) * nx + 1;
  rowEnd = i * nx;
  temp = Q;
  
  if i == N
    temp = P;
  end
  
  Qbar(rowStart:rowEnd, rowStart:rowEnd) = temp;
end

% Compute Rbar
for i = 1:N
  % R is square so we can reuse indices
  rowStart = (i - 1) * nu + 1;
  rowEnd = i * nu;
  Rbar(rowStart:rowEnd, rowStart:rowEnd) = R;
end

% Compute Optimal Control Gain
Kopt = -inv(S' * Qbar * S + Rbar) * S' * Qbar * M;
end
