clc; clear; close all;

%% system params
m = 5; % mass
k = 2; % spring coefficient
b = 0.5; % damping coefficient

Ts = 0.05; % sampling time

%% state space model
A_msd = [0 1;
  -k / m -b / m];
B_msd = [0;
  1 / m];
C_msd = [1 0];
D_msd = 0;

[Ad_msd, Bd_msd, Cd_msd, Dd_msd] = c2dm(A_msd, B_msd, C_msd, D_msd, Ts);

%% setup simulation
nx = size(Ad_msd, 1);
nu = size(Bd_msd, 2);
x0 = ones(nx, 1);
Tsim = 10;
times = 0:Ts:Tsim;
Tfid = Ts / 10;

%% deterministic

%% simulate open loop dynamics
data.x_ol = zeros(nx, length(times));
data.x_ol(:, 1) = x0;
k = 1;

for t = times(1: end - 1)
  uk = 0;
  [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t + Tfid:Tfid:t + Ts, data.x_ol(:, k));
  data.x_ol(:, k + 1) = x_ode(end, :)';
  k = k + 1;
end

%% plot open loop dynamics
fig = figure;
hold on;
plot(times, data.x_ol(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position');
% plot(times, data.x_ol(2, :), 'b', 'LineWidth', 2, 'DisplayName', 'Velocity');
title('Open loop dynamics');
xlabel('Time [s]');
ylabel('Position [m]');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/ol_det.svg');

%% set up and solve deterministic LQR problem
Q = eye(nx) * 2; % state cost
R = eye(nu) * 1; % input cost
P = eye(nx) * 10; % terminal state cost
N = length(times) - 1; % prediction horizon, excluding initial state
[Kopt, S, M, Qbar, Rbar] = solveLQR(N, Ad_msd, Bd_msd, Q, R, P);
Uopt = Kopt * x0;

%% simulate closed loop dynamics
data.x_cl = zeros(nx, length(times));
data.x_cl(:, 1) = x0;
k = 1;

for t = times(1: end - 1)
  uk = Uopt(k);
  [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t + Tfid:Tfid:t + Ts, data.x_cl(:, k));
  data.x_cl(:, k + 1) = x_ode(end, :)';
  k = k + 1;
end

%% plot closed loop dynamics
fig = figure;
hold on;
plot(times, data.x_cl(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position');
plot(times(2:end), Uopt, 'r', 'LineWidth', 2, 'DisplayName', 'Control Effort');
title('Closed loop dynamics');
xlabel('Time [s]');
ylabel('Position [m]');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/cl_det.svg');

%% stochastic with random initial state

%% set up problem with stochastic initial state
nSamples = 100;
x0_mean = 1;
x0_sd = 0.1;
x0_rv = normrnd(x0_mean, x0_sd, [nx, nSamples]);
data.x_ol_stoch = zeros(nx, length(times), nSamples);
data.x_ol_stoch(:, 1, :) = x0_rv;
Uzero = zeros(N * nu, nSamples);

%% simulate open loop dynamics
for i = 1:nSamples
  k = 1;
  for t = times(1:end - 1)
    uk = 0;
    [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t + Tfid:Tfid:t + Ts, data.x_ol_stoch(:, k, i));
    data.x_ol_stoch(:, k + 1, i) = x_ode(end, :)';
    k = k + 1;
  end
end

%% plot open loop dynamics
fig = figure;
sgtitle("Open loop dynamics with stochastic initial state (" + nSamples + " samples)");
subplot(1, 2, 1);
hold on;

for i = 1:nSamples
  plot(times, data.x_ol_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.5]);
end

plot(times, mean(data.x_ol_stoch(1, :, :), 3), 'b', 'LineWidth', 2, 'DisplayName', 'Sample Average');
xlabel('Time [s]');
ylabel('Position [m]');
legend show; legend boxoff;
grid on; grid minor;
hold off;

subplot(1, 2, 2);
hold on;
plot(times, mean(data.x_ol_stoch(1, :, :), 3), 'b', 'LineWidth', 2, 'DisplayName', 'Sample Average');
plot(times, data.x_ol(1, :), 'r', 'LineWidth', 2, 'DisplayName', 'Deterministic');
xlabel('Time [s]');
ylabel('Position [m]');
legend show; legend boxoff;
grid on; grid minor;
hold off;
fig.Position = [100 100 1000 500];
saveas(fig, 'figs/ol_stoch_init.svg');

%% set up and solve stochastic LQR problem using robust optimization
data.x_cl_stoch = zeros(nx, length(times), nSamples);
data.x_cl_stoch(:, 1, :) = x0_rv;
x0_rv_mean = mean(x0_rv, 2);
% Uopt does not depend on x0, so we can use the same Uopt
Uopt = Kopt * x0_rv_mean;

for i = 1:nSamples
  k = 1;
  for t = times(1:end - 1)
    uk = Uopt(k);
    [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t + Tfid:Tfid:t + Ts, data.x_cl_stoch(:, k, i));
    data.x_cl_stoch(:, k + 1, i) = x_ode(end, :)';
    k = k + 1;
  end
end

%% plot closed loop dynamics
fig = figure;
sgtitle("Closed loop dynamics with stochastic initial state (" + nSamples + " samples)");
subplot(1, 2, 1);
hold on;

for i = 1:nSamples
  plot(times, data.x_cl_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.5]);
end

plot(times, mean(data.x_cl_stoch(1, :, :), 3), 'b', 'LineWidth', 2, 'DisplayName', 'Sample Average');
xlabel('Time [s]');
ylabel('Position [m]');
legend show; legend boxoff;
grid on; grid minor;
hold off;

subplot(1, 2, 2);
hold on;
plot(times, mean(data.x_cl_stoch(1, :, :), 3), 'b', 'LineWidth', 2, 'DisplayName', 'Sample Average');
plot(times, data.x_cl(1, :), 'r', 'LineWidth', 2, 'DisplayName', 'Deterministic');
xlabel('Time [s]');
ylabel('Position [m]');
legend show; legend boxoff;
grid on; grid minor;
hold off;
fig.Position = [100 100 1000 500];
saveas(fig, 'figs/cl_stoch_init.svg');

%% Monte carlo estimator with variance calculation
mc_est_samples = 1000; % number of samples of each estimator for calculating the MC variance
mc_data.samples = [10, 100, 1000, 10000]; % number of samples for a single MC estimator
mc_data.var = zeros(1, length(mc_data.samples));

%% pre calculate all matrices that dont require x0
H = S' * Qbar * S + Rbar;
q_partial = S' * Qbar * M;
c_partial = M' * Qbar * M + Q;

%% run samples of monte carlo
for samples = mc_data.samples
  mc_est = zeros(mc_est_samples, 1);
  wait_bar = waitbar(0, "Running MC estimator with " + samples + " samples");
  for i = 1:mc_est_samples
    mc_x0 = normrnd(x0_mean, x0_sd, [nx, samples]);
    mc_cost = zeros(samples, 1);
    mc_x0_mean = mean(mc_x0, 2);
    Uopt = Kopt * mc_x0_mean;
    for j = 1:samples
      mc_cost(j) = Uopt' * H * Uopt + 2 * (q_partial * mc_x0(:, j))' * Uopt + mc_x0(:, j)' * c_partial * mc_x0(:, j);
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
