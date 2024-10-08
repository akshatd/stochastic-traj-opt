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
perturb_dir_samples = 100;
perturb_dir_mag = 0.1;
perturb_range = -1:0.25:6;

%% simulate open loop dynamics
[data.t, x_ode] = ode45(@(t, x) msd(t, x, 0, A_msd, B_msd), 0:Tfid:Tsim, x0);
data.x_ol = x_ode';

% plot open loop dynamics
fig = figure;
hold on;
plot(data.t, data.x_ol(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position');
plot(data.t, data.x_ol(2, :), 'g', 'LineWidth', 2, 'DisplayName', 'Velocity');
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
  plot(data.t, data.x_ol_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.2]);
end

plot(data.t, data.x_ol(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position: Deterministic', 'Color', [0 0 1, 0.5]);
plot(data.t, mean(data.x_ol_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Sample Average');
ylabel('Position [m]');
ylim(plot_ylim_l);
title("Open loop dynamics with stochastic initial state (" + rv_samples + " samples)");
xlabel('Time [s]');
legend show; legend boxoff;
grid on; grid minor;
saveas(fig, 'figs/ol_stoch_init.svg');

data.lqrsol = cell(length(TsList), 1);

for idx = 1:length(TsList)
  Ts = TsList(idx);
  fprintf("\n*** Ts: %.2f ***\n", Ts);
  assert(mod(Ts, Tfid) == 0, "Ts=" + Ts + " is not a multiple of Tfid=" + Tfid);
  
  %% deterministic LQR
  Tmult = Ts / Tfid;
  times = 0:Ts:Tsim;
  N = Tsim/Ts; % prediction horizon, excluding initial state
  % create the extended state system [x_k, u_k-1, r_k]
  x0_ext = [x0; u0; ref]; %
  [A_ext, B_ext, Q_ext, R_ext, P_ext] = extendState(c2d(msd_sys, Ts), Q, R, P);
  [Kopt, S, M, Qbar, Rbar] = solveLQR(N, A_ext, B_ext, Q_ext, R_ext, P_ext);
  Uopt = Kopt * x0_ext;
  cost_det = LQRCost(x0_ext, Uopt, Q_ext, S, M, Qbar, Rbar);
  fprintf("Deterministic LQR cost: %f\n", cost_det);
  
  %% simulate closed loop dynamics
  data.x_cl = sim_closed_loop(Tsim, Tfid, Ts, x0, u0, Uopt, A_msd, B_msd, @msd);
  
  %% plot closed loop dynamics
  fig = figure;
  hold on;
  plot(data.t, data.x_cl(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position');
  plot(data.t, data.x_cl(2, :), 'g', 'LineWidth', 2, 'DisplayName', 'Velocity');
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
  % Uopt only depends on x0, so we can use the same Kopt
  Uopt = Kopt * x0_rv_mean_ext;
  data.lqrsol{idx} = struct('x0_ext', x0_rv_mean_ext, 'Uopt', Uopt, 'Q_ext', Q_ext, 'S', S, 'M', M, 'Qbar', Qbar, 'Rbar', Rbar);
  cost_lqr = zeros(rv_samples, 1);
  x_cl_stoch = zeros(nx, Tsim/Tfid + 1, rv_samples);
  parfor i = 1:rv_samples
    cost_lqr(i) = LQRCost([x0_rv(:, i); u0; ref], Uopt, Q_ext, S, M, Qbar, Rbar);
    x_cl_stoch(:,:,i) = sim_closed_loop(Tsim, Tfid, Ts, x0_rv(:, i), u0, Uopt, A_msd, B_msd, @msd);
  end
  data.cost_lqr = cost_lqr ;
  data.x_cl_stoch = x_cl_stoch;
  
  %% plot closed loop dynamics
  fig = figure;
  hold on;
  
  for i = 1:rv_samples
    plot(data.t, data.x_cl_stoch(1, :, i), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.2]);
  end
  
  plot(data.t, data.x_cl(1, :), 'b', 'LineWidth', 2, 'DisplayName', 'Position: Deterministic', 'Color', [0 0 1, 0.5]);
  plot(data.t, mean(data.x_cl_stoch(1, :, :), 3), '--k', 'LineWidth', 2, 'DisplayName', 'Position: Sample Average');
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
  
  
  %% cost distribution
  [cost_lqr_exp, cost_lqr_var] = LQRcost_stats(x0_rv_mean_ext, x0_rv_cov_ext, Uopt, Q_ext, S, M, Qbar, Rbar);
  fprintf("Expectaion of Stochastic LQR cost\n- Analytical: %f\n- Experimental: %f\n", cost_lqr_exp, mean(data.cost_lqr));
  fprintf("Variance of Stochastic LQR cost\n- Analytical: %f\n- Experimental: %f\n", cost_lqr_var, var(data.cost_lqr));
  
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

%% correlation check for MLMC, ONLY FOR 2 LEVEL

%% low fid solution in high fid sim
% perturb the solution in the max gradient direction
perturbation = get_perturb_max_grad(data.lqrsol{1}, perturb_dir_samples, perturb_dir_mag, perturb_range);

Uopt_hf = data.lqrsol{1}.Uopt + perturbation;
% repeat low fid so it can be put into the cost fn of high fid
Uopt_lf = repmat(data.lqrsol{2}.Uopt, 10, 1) + perturbation;

% calculate costs and correlation for both high and low fidelity, for each perturbation
[cost_hf, cost_lf] = calc_costs_multifid(data.lqrsol{1}, Uopt_hf, Uopt_lf);
corr = calc_corr_multifid(x0_rv, u0, ref, data.lqrsol{1}, Uopt_hf, Uopt_lf);
%% plot costs
plot_multifid_costs(perturb_range, cost_hf, cost_lf, corr, "HF");

%% high fid solution in low fid sim
% perturb the solution in the max gradient direction
perturbation = get_perturb_max_grad(data.lqrsol{2}, perturb_dir_samples, perturb_dir_mag, perturb_range);

Uopt_hf = data.lqrsol{1}.Uopt(1:10:end, :) + perturbation; % downsample to low fid
Uopt_lf = data.lqrsol{2}.Uopt + perturbation;

% calculate costs and correlation for both high and low fidelity, for each perturbation
[cost_hf, cost_lf] = calc_costs_multifid(data.lqrsol{2}, Uopt_hf, Uopt_lf);
corr = calc_corr_multifid(x0_rv, u0, ref, data.lqrsol{2}, Uopt_hf, Uopt_lf);

%% plot costs
plot_multifid_costs(perturb_range, cost_hf, cost_lf, corr, "LF");

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

function xdot = msd(~, x, u, A, B)
xdot = A * x + B * u;
end

function cost = LQRCost(x0, u, Q, S, M, Qbar, Rbar)
cost = u'*(S'*Qbar*S + Rbar)*u + 2*x0'*M'*Qbar*S*u + x0'*(M'*Qbar*M + Q)*x0;
end

function grad = LQRgrad(x0, u, S, M, Qbar, Rbar)
% 2HU + 2q
H = S'*Qbar*S + Rbar;
q = x0'*M'*Qbar*S;
grad = 2 * H * u + 2 * q';
end


function [exp, var] = LQRcost_stats(x0_mean, x0_cov, u, Q, S, M, Qbar, Rbar)
K = u'*(S'*Qbar*S + Rbar)*u;
L = 2 * M' * Qbar * S * u;
N = M' * Qbar * M + Q;
exp = K + x0_mean' * L + x0_mean' * N * x0_mean + trace(N * x0_cov);
var = L' * x0_cov * L + 2 * trace(N * x0_cov * N * x0_cov) + 4 * (x0_mean' * N + L') * x0_cov * N * x0_mean;
end

function x = sim_closed_loop(Tsim, Tfid, Ts, x0, u0, U, A, B, dyn)
nx = size(A, 1);
x = zeros(nx, Tsim/Tfid + 1);
k = 0;
xk = x0;
ukm1 = u0;
Tmult = Ts / Tfid;
times = 0:Ts:Tsim;
for t = times(1: end - 1)
  uk = ukm1 + U(k+1);
  [~, x_ode] = ode45(@(t, x) dyn(t, x, uk, A, B), t:Tfid:t + Ts, xk);
  % skip the last x since it will be repeated in the next sim
  x(:, k*Tmult + 1:(k + 1)*Tmult) = x_ode(1:end-1,:)';
  xk = x_ode(end, :)';
  ukm1 = uk;
  k = k + 1;
end
% add back the last xk
x(:,end) = xk;
end

function perturbation = get_perturb_max_grad(lqrsol, perturb_dir_samples, perturb_dir_mag, perturb_range)
rand_grads = zeros(length(lqrsol.Uopt), perturb_dir_samples);
% generate uniform random samples in range [-1, 1]
perturb_dirs = - 1 + 2.*rand(length(lqrsol.Uopt), perturb_dir_samples);
% normalize dirs so they are unit vectors
perturb_dirs = perturb_dirs ./ sqrt(sum(perturb_dirs.^2, 1));
% scale the vector to have a magnitude of perturb_dir_mag
perturb_dirs = perturb_dirs * perturb_dir_mag;
% check the gradient in random directions
for i = 1:perturb_dir_samples
  rand_grads(:, i) = LQRgrad(lqrsol.x0_ext, lqrsol.Uopt + perturb_dirs(:, i), lqrsol.S, lqrsol.M, lqrsol.Qbar, lqrsol.Rbar);
end
% pick max gradient direction
[~, max_grad_idx] = max(sum(rand_grads.^2, 1));
perturb_dir_max = perturb_dirs(:, max_grad_idx)./ perturb_dir_mag; % length 1 so can scale later

% perturb the solution in the max gradient direction
perturbation = perturb_dir_max * perturb_range;
end

function corr = calc_corr_multifid(x0_rv, u0, ref, lqrsol, Uopt_hf, Uopt_lf)
perturbs = size(Uopt_hf, 2);
rv_samples = size(x0_rv, 2);
cost_hf_corr = zeros(perturbs, rv_samples);
cost_lf_corr = zeros(perturbs, rv_samples);
for i = 1:perturbs
  for j = 1:rv_samples
    cost_hf_corr(i, j) = LQRCost([x0_rv(:, j); u0; ref], Uopt_hf(:, i), lqrsol.Q_ext, lqrsol.S, lqrsol.M, lqrsol.Qbar, lqrsol.Rbar);
    cost_lf_corr(i, j) = LQRCost([x0_rv(:, j); u0; ref], Uopt_lf(:, i), lqrsol.Q_ext, lqrsol.S, lqrsol.M, lqrsol.Qbar, lqrsol.Rbar);
  end
end

corr = zeros(perturbs, 1);
for i = 1:perturbs
  corr_mat = corrcoef(cost_hf_corr(i, :), cost_lf_corr(i, :));
  % we only need the cross correlation, diagnonal will be 1
  corr(i) = corr_mat(1,2);
end
end

function [cost_hf, cost_lf] = calc_costs_multifid(lqrsol, Uopt_hf, Uopt_lf)
perturbs = size(Uopt_hf, 2);
cost_hf = zeros(perturbs, 1);
cost_lf = zeros(perturbs, 1);
for i = 1:perturbs
  cost_hf(i) = LQRCost(lqrsol.x0_ext, Uopt_hf(:, i), lqrsol.Q_ext, lqrsol.S, lqrsol.M, lqrsol.Qbar, lqrsol.Rbar);
  cost_lf(i) = LQRCost(lqrsol.x0_ext, Uopt_lf(:, i), lqrsol.Q_ext, lqrsol.S, lqrsol.M, lqrsol.Qbar, lqrsol.Rbar);
end
end


function plot_multifid_costs(perturb_range, cost_hf, cost_lf, corr, title_str)
fig = figure;
sgtitle("Perturbations in the " + title_str + " simulation")

subplot(2,1,1);
ax = gca;
yyaxis left % have to do this to be able to specify colors on both axes
plot(perturb_range, cost_hf, 'b', 'LineWidth', 2, 'DisplayName', 'HF solution');
ylabel('HF solution cost');
ax.YColor = 'b';
yyaxis right;
plot(perturb_range, cost_lf, 'r', 'LineWidth', 2, 'DisplayName', 'LF solution');
ylabel('LF solution cost');
ax.YColor = 'r';
title("Cost");
xlabel('Perturbation');
legend show; legend boxoff;
grid on; grid minor;

subplot(2,1,2);
plot(perturb_range, corr, 'LineWidth', 2);
title("Correlation between costs of solutions near " + title_str + " optimum")
xlabel('Perturbation');
ylabel('Correlation');
grid on; grid minor;

fig.Position(3:4) = [500, 900];
saveas(fig, "figs/perturb_comp_" + lower(title_str) + ".svg");

end