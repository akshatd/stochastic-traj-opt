clc; clear; close all force;
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

% general settings
plot_ylim_l = [-1.5, 2.5];
plot_ylim_r = [0.5, 5];

x0 = [1;0];
x0_mean = x0;
x0_cov = eye(nx)* 0.1;
ref = [2; 0]; % reference state [pos, vel]
u0 = 0; % initial control effort

Tsim = 6;
Tfid = 0.01; % simulation fidelity
% Ts = 0.05; % sampling time, MUST BE MULTIPLE OF Tfid
% Ts = 0.5; % sampling time, MUST BE MULTIPLE OF Tfid
TsList = [0.05, 0.5];

Q = eye(nx) * 2; % state cost
R = eye(nu) * 10; % input cost
% terminal state cost, needs to be 0 for the extended state since
% we dont want to penalize the state being away from the origin
P = eye(nx) * 0;

rv_samples = 1000;

% simulate open loop dynamics
[data.times, x_ode] = ode45(@(t, x) msd(t, x, 0, A_msd, B_msd), 0:Tfid:Tsim, x0);
data.x_ol = x_ode';

% plot open loop dynamics
fig = figure;
hold on;
plot(data.times, data.x_ol(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position');
plot(data.times, data.x_ol(2, :), 'g', 'LineWidth', 2, 'DisplayName', 'Velocity');
ylabel('Position [m]/Velocity [m/s]');
ylim(plot_ylim_l);
title("Open loop dynamics");
xlabel('Time [s]');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/ol_det.svg');

%% stochastic with random initial state

% set up problem with stochastic initial state
x0_rv = mvnrnd(x0_mean, x0_cov, rv_samples)';
x0_rv_mean = mean(x0_rv, 2);
x0_rv_cov = cov(x0_rv');
x0_rv_mean_ext = [x0_rv_mean; u0; ref];
x0_rv_cov_ext = blkdiag(x0_rv_cov, zeros(3,3));

% simulate open loop dynamics
data.x_ol_stoch = zeros(nx, Tsim/Tfid + 1, rv_samples);
for i = 1:rv_samples
  [~, x_ode] = ode45(@(t, x) msd(t, x, 0, A_msd, B_msd), 0:Tfid:Tsim, x0_rv(:, i));
  data.x_ol_stoch(:, :, i) = x_ode';
end

% plot open loop dynamics
fig = figure;
hold on;

for i = 1:rv_samples
  plot(data.times, data.x_ol_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.2]);
end

plot(data.times, data.x_ol(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position: Deterministic', 'Color', [0 0 1, 0.5]);
plot(data.times, mean(data.x_ol_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Sample Average');
ylabel('Position [m]');
ylim(plot_ylim_l);
title("Open loop dynamics with stochastic initial state (" + rv_samples + " samples)");
xlabel('Time [s]');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/ol_stoch_init.svg');

lqrsol = cell(length(TsList), 1);

for idx = 1:length(TsList)
  Ts = TsList(idx);
  display("Running with Ts=" + Ts + ", samples=" + rv_samples);
  assert(mod(Ts, Tfid) == 0, "Ts=" + Ts + " is not a multiple of Tfid=" + Tfid);
  
  %% deterministic LQR
  Tmult = Ts / Tfid;
  N = Tsim/Ts; % prediction horizon, excluding initial state
  % create the extended state system [x_k, u_k-1, r_k]
  x0_ext = [x0; u0; ref]; %
  [A_ext, B_ext, Q_ext, R_ext, P_ext] = extendState(c2d(msd_sys, Ts), Q, R, P);
  [Kopt, S, M, Qbar, Rbar] = solveLQR(N, A_ext, B_ext, Q_ext, R_ext, P_ext);
  Uopt = Kopt * x0_ext;
  cost_det = LQRCost(x0_ext, Uopt, Q_ext, S, M, Qbar, Rbar);
  display("Deterministic LQR cost: " + cost_det);
  
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
  ylabel('Position [m]/Velocity [m/s]');
  ylim(plot_ylim_l);
  yyaxis right;
  stairs(times(1:end-1), u0 + cumsum(Uopt), 'r', 'LineWidth', 2, 'DisplayName', 'Control Effort');
  ylabel('Control Effort [N]');
  ylim(plot_ylim_r);
  ax = gca;
  ax.YColor = 'r';
  title("(Ts:"+Ts+") Closed loop dynamics");
  xlabel('Time [s]');
  legend show; legend boxoff;
  grid on; grid minor;
  saveas(fig, "figs/"+Ts+"_cl_det.svg");
  
  %% set up and solve stochastic LQR problem using robust optimization
  data.x_cl_stoch = zeros(nx, Tsim/Tfid + 1, rv_samples);
  % Uopt only depends on x0, so we can use the same Kopt
  Uopt = Kopt * x0_rv_mean_ext;
  lqrsol{idx} = struct('x0_ext', x0_ext, 'Uopt', Uopt, 'Q_ext', Q_ext, 'S', S, 'M', M, 'Qbar', Qbar, 'Rbar', Rbar);
  data.cost_lqr = zeros(rv_samples, 1);
  for i = 1:rv_samples
    k = 0;
    xk = x0_rv(:, i);
    ukm1 = u0;
    data.cost_lqr(i) = LQRCost([x0_rv(:, i); u0; ref], Uopt, Q_ext, S, M, Qbar, Rbar);
    for t = times(1:end - 1)
      uk = ukm1 + Uopt(k+1);
      [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t:Tfid:t + Ts, xk);
      % skip the last x since it will be repeated in the next sim
      data.x_cl_stoch(:, k*Tmult + 1:(k + 1)*Tmult, i) = x_ode(1:end-1, :)';
      xk = x_ode(end, :)';
      ukm1 = uk;
      k = k + 1;
    end
    % add back the last xk
    data.x_cl_stoch(:, end, i) = xk;
  end
  
  [cost_lqr_exp, cost_lqr_var] = LQRcost_stats(x0_rv_mean_ext, x0_rv_cov_ext, Uopt, Q_ext, S, M, Qbar, Rbar);
  display("Expectaion of Stochastic LQR cost(analytical): " + cost_lqr_exp);
  display("Variance of Stochastic LQR cost(analytical): " + cost_lqr_var);
  display("Expectation of Stochastic LQR cost(MC): " + mean(data.cost_lqr));
  display("Variance of Stochastic LQR cost(MC): " + var(data.cost_lqr));
  
  %% plot closed loop dynamics
  fig = figure;
  hold on;
  
  for i = 1:rv_samples
    plot(data.times, data.x_cl_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.2]);
  end
  
  plot(data.times, data.x_cl(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position: Deterministic', 'Color', [0 0 1, 0.5]);
  plot(data.times, mean(data.x_cl_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Sample Average');
  ylabel('Position [m]');
  ylim(plot_ylim_l);
  yyaxis right;
  stairs(times(1:end-1), u0 + cumsum(Uopt), 'r', 'LineWidth', 2, 'DisplayName', 'Control Effort');
  ylabel('Control Effort [N]');
  ylim(plot_ylim_r);
  ax = gca;
  ax.YColor = 'r';
  title("(Ts:"+Ts+") Closed loop dynamics with stochastic initial state (" + rv_samples + " samples)");
  xlabel('Time [s]');
  legend show; legend boxoff;
  grid on; grid minor;
  saveas(fig, "figs/"+Ts+"_cl_stoch_init.svg");
  
  % plot cost distribution
  fig = figure;
  hold on;
  histogram(data.cost_lqr, 10, 'Normalization', 'pdf', 'FaceColor', 'b', 'EdgeColor', 'k', 'DisplayName', 'Cost Distribution');
  xline(cost_det, 'k', 'LineWidth', 2, 'DisplayName', 'Deterministic Cost');
  xline(cost_lqr_exp, 'g', 'LineWidth', 2, 'DisplayName', "Analytical Mean"+newline+"Var: "+cost_lqr_var);
  xline(mean(data.cost_lqr), '--r', 'LineWidth', 2, 'DisplayName', "Statistical Mean"+newline+"Var: "+var(data.cost_lqr));
  ylabel('Probability Density');
  title("(Ts:"+Ts+") Cost distribution of stochastic LQR ("+rv_samples+" samples)");
  xlabel('Cost');
  legend show; legend boxoff;
  grid on; grid minor;
  saveas(fig, "figs/"+Ts+"_cost_dist.svg");
end


%% perturb solutions U and check if they are still optimal
perturbation = -1:0.1:1;
for idx = 1:length(TsList)
  Ts = TsList(idx);
  sol = lqrsol{idx};
  Uopt_perturbed = sol.Uopt + perturbation;
  cost_perturbed = zeros(length(perturbation), 1);
  for i = 1:length(perturbation)
    cost_perturbed(i) = LQRCost(sol.x0_ext, Uopt_perturbed(:, i), sol.Q_ext, sol.S, sol.M, sol.Qbar, sol.Rbar);
  end
  fig = figure;
  plot(perturbation, cost_perturbed, 'b', 'LineWidth', 2, 'DisplayName', 'Cost');
  ylabel('Cost');
  title("(Ts:"+Ts+") Cost vs perturbation of U");
  xlabel('Perturbation');
  grid on; grid minor;
  saveas(fig, "figs/"+Ts+"_cost_perturb.svg");
end

%% correlation check, ONLY FOR 2 LEVEL MLMC

Uopt_hf = lqrsol{1}.Uopt + perturbation;
% repeat low fid so it can be put into the cost fn of high fid
Uopt_lf = lqrsol{2}.Uopt + perturbation;
Uopt_lf_rep = repmat(Uopt_lf, 10, 1);

% calculate costs for both high and low fidelity, for each perturbation
cost_hf = zeros(length(perturbation), 1);
cost_lf = zeros(length(perturbation), 1);
for i = 1:length(perturbation)
  cost_hf(i) = LQRCost(lqrsol{1}.x0_ext, Uopt_hf(:, i), lqrsol{1}.Q_ext, lqrsol{1}.S, lqrsol{1}.M, lqrsol{1}.Qbar, lqrsol{1}.Rbar);
  cost_lf(i) = LQRCost(lqrsol{1}.x0_ext, Uopt_lf_rep(:, i), lqrsol{1}.Q_ext, lqrsol{1}.S, lqrsol{1}.M, lqrsol{1}.Qbar, lqrsol{1}.Rbar);
end
% plot both costs
fig = figure;
ax = gca;
yyaxis left % have to do this to be able to specify colors on both axes
plot(perturbation, cost_hf, 'b', 'LineWidth', 2, 'DisplayName', 'Cost HF');
ylabel('High fidelity solution');
ax.YColor = 'b';
yyaxis right;
plot(perturbation, cost_lf, 'r', 'LineWidth', 2, 'DisplayName', 'Cost LF');
ylabel('Low fidelity solution');
ax.YColor = 'r';
title("HF simulation cost vs HF and LF solutions with perturbation");
xlabel('Perturbation');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, "figs/cost_perturb_comp_hf.svg");

%% Monte carlo estimator with variance calculation
mc_data.samples = [10, 100, 1000, 10000]; % number of samples for a single MC estimator
mc_data.var = zeros(1, length(mc_data.samples));

% go through MC estimators with different sample sizes
for samples = mc_data.samples
  mc_est = zeros(rv_samples, 1);
  wait_bar = waitbar(0, "Running MC estimator with " + samples + " samples");
  % run each MC estimator multiple times to get a vaiance estimate
  for i = 1:rv_samples
    x0_mc = mvnrnd(x0_mean, x0_cov, samples)'; % take samples for a single estimator
    x0_mc_mean = mean(x0_mc, 2);
    x0_mc_mean_ext = [x0_mc_mean; u0; ref];
    Uopt = Kopt * x0_mc_mean_ext; % control effort for the mean initial state
    mc_x = zeros(nx, Tsim/Tfid + 1, samples);
    % simulate trajectories for each sample
    for j = 1:samples
      k = 0;
      xk = x0_mc(:, j);
      ukm1 = u0;
      for t = times(1:end - 1)
        uk = ukm1 + Uopt(k+1);
        [~, x_ode] = ode45(@(t, x) msd(t, x, uk, A_msd, B_msd), t:Tfid:t + Ts, xk);
        % skip the last x since it will be repeated in the next sim
        mc_x(:, k*Tmult + 1:(k + 1)*Tmult, j) = x_ode(1:end-1, :)';
        xk = x_ode(end, :)';
        ukm1 = uk;
        k = k + 1;
      end
      % add back the last xk
      mc_x(:, end, j) = xk;
    end
    % calculate the MSE of the trajectory
    mc_cost = squeeze(sum((mc_x - ref).^2, [1,2]));
    
    waitbar(i / rv_samples, wait_bar);
    mc_est(i) = mean(mc_cost);
  end
  close(wait_bar);
  mc_data.var(mc_data.samples == samples) = var(mc_est);
end

%% plot MC variance
fig = figure();
loglog(mc_data.samples, mc_data.var, 'b', 'LineWidth', 2, 'DisplayName', 'MC Variance');
hold on;
title("(Ts:"+Ts+") MC estimator variance of mean sq error");
xlabel('Number of samples');
ylabel('Variance');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, "figs/"+Ts+"_mc_variance.svg");
hold off;

function xdot = msd(t, x, u, A, B)
xdot = A * x + B * u;
end

function cost = LQRCost(x0, u, Q, S, M, Qbar, Rbar)
cost = u'*(S'*Qbar*S + Rbar)*u + 2*x0'*M'*Qbar*S*u + x0'*(M'*Qbar*M + Q)*x0;
end

function [exp, var] = LQRcost_stats(x0_mean, x0_cov, u, Q, S, M, Qbar, Rbar)
K = u'*(S'*Qbar*S + Rbar)*u;
L = 2 * M' * Qbar * S * u;
N = M' * Qbar * M + Q;
exp = K + x0_mean' * L + x0_mean' * N * x0_mean + trace(N * x0_cov);
var = L' * x0_cov * L + 2 * trace(N * x0_cov * N * x0_cov) + 4 * (x0_mean' * N + L') * x0_cov * N * x0_mean;
end