clc; clear; close all;

%% state space model

% system params
m = 5; % mass
k = 2; % spring coefficient
b = 0.5; % damping coefficient

A_msd = [0 1;
  -k / m -b / m];
B_msd = [0;
  1 / m];
% C_msd = [1 0]; % only observe position
C_msd = [1 0;
  0 1]; % observe both position and velocity
D_msd = 0;

msd_sys = ss(A_msd, B_msd, C_msd, D_msd);

% setup simulation
nx = size(A_msd, 1);
nu = size(B_msd, 2);
x0 = [1;0];
Tsim = 10;
Tfid = 0.01; % simulation fidelity

% simulate open loop dynamics
[data.times, x_ode] = ode45(@(t, x) msd(t, x, 0, A_msd, B_msd), 0:Tfid:Tsim, x0);
data.x_ol = x_ode';

% plot open loop dynamics
fig = figure;
hold on;
plot(data.times, data.x_ol(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position');
plot(data.times, data.x_ol(2, :), 'g', 'LineWidth', 2, 'DisplayName', 'Velocity');
title('Open loop dynamics');
xlabel('Time [s]');
ylabel('Position [m]/Velocity [m/s]');
ylim([-3 3]);
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/ol_det.svg');

%% deterministic LQR
Ts = 0.05; % sampling time, MUST BE MULTIPLE OF Tfid
assert(mod(Ts, Tfid) == 0, "Ts=" + Ts + " is not a multiple of Tfid=" + Tfid);
Tmult = Ts / Tfid;
%% set up and solve deterministic LQR problem
Q = eye(nx) * 2; % state cost
R = eye(nu) * 100; % input cost
% terminal state cost, needs to be 0 for the extended state since
% we dont want to penalize the state being away from the origin
P = eye(nx) * 0;
N = Tsim/Ts; % prediction horizon, excluding initial state
ref = [-2; 0]; % reference state [pos, vel]
u0 = 0; % initial control effort
% create the extended state system [x_k, u_k-1, r_k]
x0_ext = [x0; u0; ref]; %
[A_ext, B_ext, Q_ext, R_ext, P_ext] = extendState(c2d(msd_sys, Ts), Q, R, P);
Kopt = solveLQR(N, A_ext, B_ext, Q_ext, R_ext, P_ext);
Uopt = Kopt * x0_ext;

%% simulate closed loop dynamics
data.x_cl = zeros(nx, Tsim/Tfid + 1);
k = 0;
xk = x0;
ukm1 = u0;
times = 0:Ts:Tsim;
for t = times(1: end - 1)
  uk = ukm1 + Uopt(k+1);
  [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t:Tfid:t + Ts, xk);
  % skip the last x since it will be repeated in the next sim
  data.x_cl(:, k*Tmult + 1:(k + 1)*Tmult) = x_ode(1:end-1,:)';
  xk = x_ode(end, :)';
  ukm1 = uk;
  k = k + 1;
end
% add back the last xk
data.x_cl(:,end) = xk;

%% plot closed loop dynamics
fig = figure;
hold on;
plot(data.times, data.x_cl(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position');
plot(data.times, data.x_cl(2, :), 'g', 'LineWidth', 2, 'DisplayName', 'Velocity');
stairs(times(1:end-1), u0 + cumsum(Uopt), 'r', 'LineWidth', 2, 'DisplayName', 'Control Effort');
title('Closed loop dynamics');
xlabel('Time [s]');
ylabel('Position [m]/Velocity [m/s]');
ylim([-3 3]);
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/cl_det.svg');

%% stochastic with random initial state

%% set up problem with stochastic initial state
nSamples = 100;
x0_mean = x0;
x0_sd = eye(nx)* 0.1;
x0_rv = mvnrnd(x0_mean, x0_sd, nSamples)';
data.x_ol_stoch = zeros(nx, Tsim/Tfid + 1, nSamples);

%% simulate open loop dynamics
for i = 1:nSamples
  [~, x_ode] = ode45(@(t, x) msd(t, x, 0, A_msd, B_msd), 0:Tfid:Tsim, x0_rv(:, i));
  data.x_ol_stoch(:, :, i) = x_ode';
end

%% plot open loop dynamics
fig = figure;
sgtitle("Open loop dynamics with stochastic initial state (" + nSamples + " samples)");
subplot(1, 2, 1);
hold on;

for i = 1:nSamples
  plot(data.times, data.x_ol_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.5]);
end

plot(data.times, mean(data.x_ol_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Position: Sample Average');
xlabel('Time [s]');
ylabel('Position [m]');
ylim([-3 3]);
legend show; legend boxoff;
grid on; grid minor;
hold off;

subplot(1, 2, 2);
hold on;
plot(data.times, data.x_ol(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position: Deterministic');
plot(data.times, mean(data.x_ol_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Position: Sample Average');
xlabel('Time [s]');
ylabel('Position [m]');
ylim([-3 3]);
legend show; legend boxoff;
grid on; grid minor;
hold off;
fig.Position = [100 100 1000 500];
saveas(fig, 'figs/ol_stoch_init.svg');

%% set up and solve stochastic LQR problem using robust optimization
data.x_cl_stoch = zeros(nx, Tsim/Tfid + 1, nSamples);
x0_rv_mean = mean(x0_rv, 2);
% Uopt does not depend on x0, so we can use the same Uopt
Uopt = Kopt * x0_rv_mean;

for i = 1:nSamples
  k = 0;
  xk = x0_rv(:, i);
  for t = times(1:end - 1)
    uk = Uopt(k+1);
    [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t:Tfid:t + Ts, xk);
    % skip the last x since it will be repeated in the next sim
    data.x_cl_stoch(:, k*Tmult + 1:(k + 1)*Tmult, i) = x_ode(1:end-1, :)';
    xk = x_ode(end, :)';
    k = k + 1;
  end
  % add back the last xk
  data.x_cl_stoch(:, end, i) = xk;
end

%% plot closed loop dynamics
fig = figure;
sgtitle("Closed loop dynamics with stochastic initial state (" + nSamples + " samples)");
subplot(1, 2, 1);
hold on;

for i = 1:nSamples
  plot(data.times, data.x_cl_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.5]);
end

plot(data.times, mean(data.x_cl_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Sample Average');
stairs(times(1:end-1), Uopt, 'r', 'LineWidth', 2, 'DisplayName', 'Control Effort');
xlabel('Time [s]');
ylabel('Position [m]');
ylim([-3 3]);
legend show; legend boxoff;
grid on; grid minor;
hold off;

subplot(1, 2, 2);
hold on;
plot(data.times, data.x_cl(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position: Deterministic');
plot(data.times, mean(data.x_cl_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Sample Average');
xlabel('Time [s]');
ylabel('Position [m]');
ylim([-3 3]);
legend show; legend boxoff;
grid on; grid minor;
hold off;
fig.Position = [100 100 1000 500];
saveas(fig, 'figs/cl_stoch_init.svg');

%% Monte carlo estimator with variance calculation
mc_est_samples = 1000; % number of samples of each estimator for calculating the MC variance
mc_data.samples = [10, 100, 1000, 10000]; % number of samples for a single MC estimator
mc_data.var = zeros(1, length(mc_data.samples));

% go through MC estimators with different sample sizes
for samples = mc_data.samples
  mc_est = zeros(mc_est_samples, 1);
  wait_bar = waitbar(0, "Running MC estimator with " + samples + " samples");
  % run each MC estimator multiple times to get a vaiance estimate
  for i = 1:mc_est_samples
    x0_mc = normrnd(x0_mean, x0_sd, [nx, samples]);
    mc_cost = zeros(samples, 1);
    x0_mc_mean = mean(x0_mc, 2);
    Uopt = Kopt * x0_mc_mean;
    mc_x = zeros(nx, Tsim/Tfid + 1, samples);
    % for j = 1:samples
    %   mc_cost(j) = Uopt' * H * Uopt + 2 * (q_partial * x0_mc(:, j))' * Uopt + x0_mc(:, j)' * c_partial * x0_mc(:, j);
    % end
    % calculate the cost for each sample
    for j = 1:samples
      k = 0;
      xk = x0_mc(:, j);
      for t = times(1:end - 1)
        uk = Uopt(k+1);
        [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t:Tfid:t + Ts, xk);
        % skip the last x since it will be repeated in the next sim
        mc_x(:, k*Tmult + 1:(k + 1)*Tmult, j) = x_ode(1:end-1, :)';
        xk = x_ode(end, :)';
        k = k + 1;
      end
      % add back the last xk
      mc_x(:, end, j) = xk;
      % alculate the total cost
      % for t=1:Tsim/Tfid
    end
    
    waitbar(i / mc_est_samples, wait_bar);
    mc_est(i) = mean(mc_cost) / N;
  end
  close(wait_bar);
  mc_data.var(mc_data.samples == samples) = var(mc_est);
end

%% plot MC variance
fig = figure();
loglog(mc_data.samples, mc_data.var, 'b', 'LineWidth', 2, 'DisplayName', 'MC Variance');
hold on;
title('Monte Carlo estimator variance of normalized cost function');
xlabel('Number of samples');
ylabel('Variance');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/mc_variance.svg');
hold off;

function xdot = msd(t, x, u, A, B)
xdot = A * x + B * u;
end
