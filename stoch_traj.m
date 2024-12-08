clc; clear; close all force;
%% A. state space model

% system params
m = 5; % mass
k = 2; % spring coefficient
c = 0.5; % damping coefficient

A_msd = [0 1;
  -k / m -c / m];
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

rv_samples = 10000;
perturb_dir_samples = 1000;
perturb_dir_mag = 0.1;
perturb_range = -1:0.25:6;

%% B. open loop dynamics
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
legend show;
grid on;
saveas(fig, 'figs/ol_det.svg');

%% C. stochastic dynamics with random initial state

% set up problem with stochastic initial state
x0_rv = mvnrnd(x0_mean, x0_cov, rv_samples)';
x0_rv_mean = mean(x0_rv, 2);
x0_rv_cov = cov(x0_rv');
x0_rv_mean_ext = [x0_rv_mean; u0; ref];
x0_rv_cov_ext = blkdiag(x0_rv_cov, zeros(3,3));

%% C.1 open loop dynamics
data.x_ol_stoch = zeros(nx, Tsim/Tfid + 1, rv_samples);
for i = 1:rv_samples
  [~, x_ode] = ode45(@(t, x) msd(t, x, 0, A_msd, B_msd), 0:Tfid:Tsim, x0_rv(:, i));
  data.x_ol_stoch(:, :, i) = x_ode';
end

%% C.2 plot open loop dynamics
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
legend show;
grid on;
saveas(fig, 'figs/ol_stoch_init.svg');

%% D closed loop dynamics with LQR

data.lqrsol = cell(length(TsList), 1);

for idx = 1:length(TsList)
  Ts = TsList(idx);
  fprintf("\n*** Ts: %.2f ***\n", Ts);
  assert(mod(Ts, Tfid) == 0, "Ts=" + Ts + " is not a multiple of Tfid=" + Tfid);
  
  %% D.1 deterministic LQR
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
  legend show;
  grid on;
  saveas(fig, "figs/"+Ts+"_cl_det.svg");
  
  %% D.2 stochastic LQR using robust SAA
  % Uopt only depends on x0, so we can use the same Kopt
  Uopt = Kopt * x0_rv_mean_ext;
  data.lqrsol{idx} = struct('x0_ext', x0_rv_mean_ext, 'Uopt', Uopt, 'Q_ext', Q_ext, 'S', S, 'M', M, 'Qbar', Qbar, 'Rbar', Rbar, 'times', times);
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
  legend show;
  grid on;
  saveas(fig, "figs/"+Ts+"_cl_stoch_init.svg");
  
  
  %% D.3 cost distribution
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
  legend show;
  grid on;
  saveas(fig, "figs/"+Ts+"_cost_dist.svg");
end

% *** E. correlation check for MLMC, ONLY FOR 2 LEVEL TODO: extend to work for any number of levels
% E.1 correlation with perturbation in direction of maximum ascent

%% E.1.1 perturbation in HF objective fn
perturbation = get_perturb_max_grad(data.lqrsol{1}, perturb_dir_samples, perturb_dir_mag, perturb_range);
Uopt_hf = data.lqrsol{1}.Uopt + perturbation; % perturb hf sol

%% E.1.1.1 restriction
Uopt_lf = Uopt_hf(1:10:end, :); % pick every 10th elem
% Uopt_lf_vis = repelem(Uopt_lf, 10, 1); % uncomment to visualize
% vis_sols(Uopt_hf, Uopt_lf_vis, data.lqrsol{1}.times, "HF obj", perturb_range, "Perturbation");

% calculate costs and correlation for HF/LF sols for each perturbation
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
[corr_st, corr_an] = calc_corr(cost_hf, cost_lf, x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
%% plot costs
title_str = ["$J_h$ and $J_l$", "$J_h$ (restriction)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hlr})$"];
plot_multifid_costs(perturb_range, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Perturbation");

%% E.1.1.2 averaging
Uopt_lf = downsample_avg(Uopt_hf, 10); % get HF in LF by avg every 10 values of HF
% vis_sols(Uopt_hf, Uopt_lf, data.lqrsol{1}.times, "HF obj", perturb_range, "Perturbation");
Uopt_lf = Uopt_lf(1:10:end, :); % pick every 10th elem
% calculate costs and correlation for HF/LF sols for each perturbation
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
[corr_st, corr_an] = calc_corr(cost_hf, cost_lf, x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
% plot costs
title_str = ["$J_h$ and $J_l$", "$J_h (averaging)$"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla})$"];
plot_multifid_costs(perturb_range, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Perturbation");


%% E.1.2 perturbation in LF objective fn
perturbation = get_perturb_max_grad(data.lqrsol{2}, perturb_dir_samples, perturb_dir_mag, perturb_range);
Uopt_lf = data.lqrsol{2}.Uopt + perturbation;
% vis_sols(Uopt_hf, Uopt_lf, data.lqrsol{2}.times, "LF obj (restriction)", perturb_range, "Perturbation");
Uopt_hf = repelem(Uopt_lf, 10, 1); % repeat the LF sol to match HF sol
% calculate costs and correlation for both high and low fidelity, for each perturbation
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
[corr_st, corr_an] = calc_corr(cost_hf, cost_lf, x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
% plot costs
title_str = ["$J_h$ and $J_l$", "$J_l$"];
costs_str = ["$J_h(u_lh)$", "$J_l(u_l)$"];
plot_multifid_costs(perturb_range, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Perturbation");

%% *** E.2 correlation at different points in a numerical optimizer
%% E.2.1 num opt with HF objective fn
Uopt_hf = data.lqrsol{1}.Uopt;
% u0_num = zeros(size(Uopt_hf));
u0_num = repelem(data.lqrsol{2}.Uopt, 10, 1); % warm start
global num_opt_iters num_opt_data;
num_opt_iters = 0;
num_opt_data = zeros(size(Uopt_hf,1), 100);
fun = @(u) LQRCostwithGrad(data.lqrsol{1}.x0_ext, u, data.lqrsol{1}.Q_ext, data.lqrsol{1}.S, data.lqrsol{1}.M, data.lqrsol{1}.Qbar, data.lqrsol{1}.Rbar);
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, 'OutputFcn', @outfun);
Uopt_num = fminunc(fun, u0_num, options);
uopt_diff = sqrt(sum((Uopt_hf - Uopt_num).^2));

U_hf = num_opt_data(:, 1:num_opt_iters);

% this is just to get the costs for the averaged HF solution in LF sim as the MLMC estimator would
%% E.2.1.1 correlation with restriction
U_lf = U_hf(1:10:end, :);
% U_lf = repelem(U_lf, 10, 1); % uncomment to visualize
% vis_sols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:num_opt_iters, "Iteration");
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
[corr_st, corr_an] = calc_corr(cost_hf, cost_lf, x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);

title_str = ["$J_h$ and $J_l$", "$J_h$ (restriction)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hlr})$"];
plot_multifid_costs(1:num_opt_iters, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Iteration");
% Get correlation for all iterations
corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
%% plot correlation
plot_corr_2d(corr, "Correlation between costs in $J_h(u_h)$ and $J_l(u_{hlr})$", costs_str, "num_iters_res");

%% E.2.1.2 correlation with averaging
U_lf = downsample_avg(U_hf, 10); % get LF by avg every 10 values of HF
% vis_sols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:num_opt_iters, "Iteration");
U_lf = U_lf(1:10:end, :);
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
[corr_st, corr_an] = calc_corr(cost_hf, cost_lf, x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);

title_str = ["$J_h$ and $J_l$", "$J_h$ (averaging)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla})$"];
plot_multifid_costs(1:num_opt_iters, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Iteration");
corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
%% plot correlation
plot_corr_2d(corr, "Correlation between costs in $J_h(u_h)$ and $J_l(u_{hla})$", costs_str, "num_iters_avg");
% save all data
save("artefacts/data.mat");

%% num opt with MLMC estimator
% safe case, alpha of 0 to eliminate LF evaluations
alpha = 0.5;
num_opt_iters = 0;
num_opt_data = zeros(size(Uopt_hf,1), 100);
fun = @(u) calc_mlmc_est(data.lqrsol{1}, data.lqrsol{2}, u, 100, 100, u0, ref, alpha, x0_rv);
% fun = @(u) calc_mlmc_est(data.lqrsol{1}, data.lqrsol{2}, u, 100, 1, u0, ref, alpha, x0_mean, x0_cov);
options = optimoptions('fminunc', 'Display', 'iter', 'OptimalityTolerance', 1e-3, 'OutputFcn', @outfun);
Uopt_num = fminunc(fun, u0_num, options);
U_hf = num_opt_data(:, 1:num_opt_iters);
U_lf = downsample_avg(U_hf, 10);
vis_sols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:num_opt_iters, "Iteration");

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
legend show;
grid on;
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

function [f, g] = LQRCostwithGrad(x0, u, Q, S, M, Qbar, Rbar)
f = LQRCost(x0, u, Q, S, M, Qbar, Rbar);
g = LQRgrad(x0, u, S, M, Qbar, Rbar);
end

function stop = outfun(x, ~, state)
% this is for capturing intermediate values of the optimization
stop = false;
global num_opt_iters num_opt_data;
if isequal(state, 'iter')
  num_opt_iters = num_opt_iters + 1;
  num_opt_data(:, num_opt_iters) = x;
end
end

function [exp, var] = LQRcost_stats(x0_mean, x0_cov, u, Q, S, M, Qbar, Rbar)
K = u'*(S'*Qbar*S + Rbar)*u;
L = 2 * M' * Qbar * S * u;
N = M' * Qbar * M + Q;
exp = K + x0_mean' * L + x0_mean' * N * x0_mean + trace(N * x0_cov);
var = L' * x0_cov * L + 2 * trace(N * x0_cov * N * x0_cov) + 4 * (x0_mean' * N + L') * x0_cov * N * x0_mean;
end

function cov_lqr = LQRcost_cov(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, Uopt_hf, Uopt_lf)
% K1 = Uopt_hf' * (lqrsol_hf.S' * lqrsol_hf.Qbar * lqrsol_hf.S + lqrsol_hf.Rbar) * Uopt_hf;
L1 = 2 * lqrsol_hf.M' * lqrsol_hf.Qbar * lqrsol_hf.S * Uopt_hf;
N1 = lqrsol_hf.M' * lqrsol_hf.Qbar * lqrsol_hf.M + lqrsol_hf.Q_ext;
% K2 = Uopt_lf' * (lqrsol_lf.S' * lqrsol_lf.Qbar * lqrsol_lf.S + lqrsol_lf.Rbar) * Uopt_lf;
L2 = 2 * lqrsol_lf.M' * lqrsol_lf.Qbar * lqrsol_lf.S * Uopt_lf;
N2 = lqrsol_lf.M' * lqrsol_lf.Qbar * lqrsol_lf.M + lqrsol_lf.Q_ext;
x0_rv_ext = [x0_rv; repmat(u0, 1, size(x0_rv, 2)); repmat(ref, 1, size(x0_rv, 2))];
x_mean = mean(x0_rv_ext, 2);
x_cov = cov(x0_rv_ext');
cov_lqr = L1'*x_cov*L2 + 2*x_mean'*N2*x_cov*L1 + 2*x_mean'*N1*x_cov*L2 + 2*trace(N1*x_cov*N2*x_cov) + 4*x_mean'*N1*x_cov*N2*x_mean;
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

function U_lf = downsample_avg(U_hf, avg_window)
U_hf_len = size(U_hf, 1);
U_lf = zeros(U_hf_len/avg_window, size(U_hf, 2));
skips = U_hf_len/avg_window;
for i = 1:skips
  row_start = (i-1) * avg_window + 1;
  row_end = i*avg_window;
  U_lf(i, :) = mean(U_hf(row_start:row_end, :), 1);
end
U_lf = repelem(U_lf, avg_window, 1);
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

function [corr_st, corr_an] = calc_corr(cost_hf, cost_lf, x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, U_hf, U_lf)
perturbs = size(cost_hf, 1);
corr_st = zeros(perturbs, 1);
corr_an = zeros(perturbs, 1);
x0_rv_ext = [x0_rv; repmat(u0, 1, size(x0_rv, 2)); repmat(ref, 1, size(x0_rv, 2))];
x_mean = mean(x0_rv_ext, 2);
x_cov = cov(x0_rv_ext');
for i = 1:perturbs
  corr_mat = corrcoef(cost_hf(i, :), cost_lf(i, :));
  corr_st(i) = corr_mat(1,2); % we only need the cross correlation, diagnonal will be 1
  cov_an = LQRcost_cov(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, U_hf(:, i), U_lf(:, i)); % TODO make indices consistent
  [~, var_J1] = LQRcost_stats(x_mean, x_cov, U_hf(:, i), lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
  [~, var_J2] = LQRcost_stats(x_mean, x_cov, U_lf(:, i), lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
  corr_an(i) = cov_an / sqrt(var_J1 * var_J2);
end
end

function corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, Uopt_hf, Uopt_lf)
iters = size(Uopt_hf, 2);
rv_samples = size(x0_rv, 2);
cost_hf_corr = zeros(rv_samples, iters);
cost_lf_corr = zeros(rv_samples, iters);
for i = 1:iters
  for j = 1:rv_samples
    cost_hf_corr(j, i) = LQRCost([x0_rv(:, j); u0; ref], Uopt_hf(:, i), lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
    cost_lf_corr(j, i) = LQRCost([x0_rv(:, j); u0; ref], Uopt_lf(:, i), lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
  end
end

corr = zeros(iters, iters);
for i = 1:iters % HF iterations
  for j = 1:iters % LF iterations
    corr_mat = corrcoef(cost_hf_corr(:, i), cost_lf_corr(:, j));
    % we only need the cross correlation, diagnonal will be 1
    corr(i, j) = corr_mat(1,2); % rowas are HF, cols are LF
  end
end
end

function [cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, Uopt_hf, Uopt_lf)
perturbs = size(Uopt_hf, 2);
rv_samples = size(x0_rv, 2);
cost_hf = zeros(perturbs, rv_samples);
cost_lf = zeros(perturbs, rv_samples);
for i = 1:perturbs
  for j = 1:rv_samples
    cost_hf(i,j) = LQRCost([x0_rv(:, j); u0; ref], Uopt_hf(:, i), lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
    cost_lf(i,j) = LQRCost([x0_rv(:, j); u0; ref], Uopt_lf(:, i), lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
  end
end
end

function plot_multifid_costs(perturb_range, cost_hf, cost_lf, corr_st, corr_an, title_str, costs_str, xaxis_str)
fig = figure;
sgtitle("Mean cost in " + title_str(1) + " vs " + xaxis_str + " in " + title_str(2), 'Interpreter', 'latex')

subplot(2,1,1);
title("Mean cost");
xlabel(xaxis_str);
ax = gca;
yyaxis left % have to do this to be able to specify colors on both axes
semilogy(perturb_range, cost_hf, 'b', 'LineWidth', 2, 'DisplayName', costs_str(1));
ylabel(costs_str(1) + " cost", "Interpreter", "latex", "FontSize", 14);
ax.YColor = 'b';
yyaxis right;
semilogy(perturb_range, cost_lf, 'r', 'LineWidth', 2, 'DisplayName', costs_str(2));
ylabel(costs_str(2) + " cost", "Interpreter", "latex", "FontSize", 14);
ax.YColor = 'r';
legend('Interpreter', 'latex', 'Location', 'best', "FontSize", 12);
grid on;

subplot(2,1,2);
plot(perturb_range, corr_an, 'LineWidth', 2, 'DisplayName', 'Analytical');
hold on;
plot(perturb_range, corr_st, '--', 'LineWidth', 2, 'DisplayName', 'Statistical');
title("Correlation");
xlabel(xaxis_str);
ylabel('Pearson Correlation coefficient');
legend show;
grid on;

fig.Position(3:4) = [500, 900];
saveas(fig, "figs/perturb_comp_" + lower(title_str(1)) + "_" + title_str(2) + ".svg");

end

function plot_corr_2d(corr, title_str, sols_str, save_str)
fig = figure;
sgtitle(title_str, 'Interpreter', 'latex');

subplot(1,2,1);
imagesc(corr);
set(gca, 'YDir', 'normal');
title("Heatmap");
ylabel(sols_str(1) + " at iteration", "Interpreter", "latex", "FontSize", 14);
xlabel(sols_str(2) + " at iteration", "Interpreter", "latex", "FontSize", 14);
colorbar;
grid on;

subplot(1,2,2);
surf(corr);
title("3D plot");
ylabel(sols_str(1) + " at iteration", "Interpreter", "latex", "FontSize", 14);
xlabel(sols_str(2) + " at iteration", "Interpreter", "latex", "FontSize", 14);
zlabel("Correlation");
view(3);
grid on;

fig.Position(3:4) = [1500, 600];
saveas(fig, "figs/" + save_str + ".svg");
end

function vis_sols(uopt_hf, uopt_lf, times, title_str, zaxis, zaxis_str)
fig = uifigure("Name", title_str);
fig.Position(3:4) = [500, 800];
grd = uigridlayout(fig, [6, 1], "RowHeight", {'1x', '10x', '1x', '10x', '1x', 50});

title = uilabel(grd,'Text','Control output U');
title.Layout.Row = 1;
title.HorizontalAlignment = 'center';
title.FontSize = 16;

ax = uiaxes(grd);
ax.Layout.Row = 2;
hold(ax, 'on');
% TODO: we are assuming u0 to be 0, might not be true in general so might need to add u0 to cumsum
plot_hf = stairs(ax, times(1:end-1), cumsum(uopt_hf(:, 1)), 'b', 'LineWidth', 2, 'DisplayName', 'HF solution');
plot_lf = stairs(ax, times(1:end-1), cumsum(uopt_lf(:, 1)), 'r', 'LineWidth', 2, 'DisplayName', 'LF solution');
xlabel(ax, 'Time [s]');
ylabel(ax, 'Control Effort [N]');
grid(ax, 'on')
legend(ax, 'Location', 'best', 'Box', 'off');
hold(ax, 'off');

title = uilabel(grd,'Text','Change in control output U');
title.Layout.Row = 3;
title.HorizontalAlignment = 'center';
title.FontSize = 16;

ax = uiaxes(grd);
ax.Layout.Row = 4;
hold(ax, 'on');
plot_hf_raw = stairs(ax, times(1:end-1), uopt_hf(:, 1), 'b', 'LineWidth', 2, 'DisplayName', 'HF solution');
plot_lf_raw = stairs(ax, times(1:end-1), uopt_lf(:, 1), 'r', 'LineWidth', 2, 'DisplayName', 'LF solution');
xlabel(ax, 'Time [s]');
ylabel(ax, 'Control Effort [N]');
grid(ax, 'on')
legend(ax, 'Location', 'best', 'Box', 'off');
hold(ax, 'off');

title = uilabel(grd,'Text',zaxis_str);
title.Layout.Row = 5;
title.HorizontalAlignment = 'center';
title.FontSize = 16;

sld = uislider(grd, 'slider', 'limits', [zaxis(1), zaxis(end)]);
sld.Layout.Row = 6;
sld.ValueChangedFcn = @(sld, event) update_plot(sld, plot_hf, plot_lf, plot_hf_raw, plot_lf_raw, uopt_hf, uopt_lf);
end

function update_plot(sld, plot_hf, plot_lf, plot_hf_raw, plot_lf_raw, uopt_hf, uopt_lf)
% handle floats in the zaxis and map them back to the data range
idx = round(interp1(sld.Limits, [1, size(uopt_hf, 2)], sld.Value));

hold(plot_hf.Parent, 'on');
plot_hf.YData = cumsum(uopt_hf(:, idx));
plot_lf.YData = cumsum(uopt_lf(:, idx));
hold(plot_hf.Parent, 'off');

hold(plot_hf_raw.Parent, 'on');
plot_hf_raw.YData = uopt_hf(:, idx);
plot_lf_raw.YData = uopt_lf(:, idx);
hold(plot_hf_raw.Parent, 'off');
end

function cost = calc_mlmc_est(lqrsol_hf, lqrsol_lf, u_hf, n_hf, n_lf, u0, ref, alpha, x0_rv)
cost_hf = 0;
for i=1:n_hf
  cost_hf = cost_hf + LQRCost([x0_rv(:, i); u0; ref], u_hf, lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
end
cost_hf = cost_hf / n_hf;

u_lf = downsample_avg(u_hf, 10);
u_lf = u_lf(1:10:end);
cost_lf = 0;
% reuse same samples for LF
for i=1:n_hf
  cost_lf = cost_lf + LQRCost([x0_rv(:, i); u0; ref], u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
end
cost_lf = cost_lf / n_hf;

exp_lf = 0;
for i=n_hf:n_hf+n_lf
  exp_lf = exp_lf + LQRCost([x0_rv(:, i); u0; ref], u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
end
exp_lf = exp_lf/n_lf;

cost = cost_hf + alpha * (cost_lf - exp_lf);

end

% function cost = calc_mlmc_est(lqrsol_hf, lqrsol_lf, u_hf, n_hf, n_lf, u0, ref, alpha, x0_mean, x0_cov)
% x0_rv_hf = mvnrnd(x0_mean, x0_cov, n_hf)';
% cost_hf = 0;
% for i=1:n_hf
%   cost_hf = cost_hf + LQRCost([x0_rv_hf(:, i); u0; ref], u_hf, lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
% end
% cost_hf = cost_hf / n_hf;
%
% u_lf = downsample_avg(u_hf, 10);
% u_lf = u_lf(1:10:end);
% cost_lf = 0;
% % reuse same samples for LF
% for i=1:n_hf
%   cost_lf = cost_lf + LQRCost([x0_rv_hf(:, i); u0; ref], u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
% end
% cost_lf = cost_lf / n_hf;
%
% % Expectation of LF
% x0_rv_lf_exp = mvnrnd(x0_mean, x0_cov, n_lf)';
% exp_lf = 0;
% for i=1:n_lf
%   exp_lf = exp_lf + LQRCost([x0_rv_lf_exp(:, i); u0; ref], u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
% end
% exp_lf = exp_lf / n_lf;
%
% % cost = cost_hf + alpha * (cost_lf - exp_lf);
% cost = cost_hf;
%
% end