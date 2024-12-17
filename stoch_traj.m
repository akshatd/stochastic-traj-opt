% TODOS:
% - we are assuming u0 to be 0, might not be true in general so might need to add u0 to cumsum
% - extend MLMC to work for any number of levels
% - make indices consistent when storing U in row vs column (make everything row based)
% - split out all the stats related functions into a separate file
% - calculate correlation in cv estimator as it goes along instead of all at once
% - calculate 2d correlation analytically
% - make U and u consistent in naming
% - make usage of x_ext vs x, u and ref consistent, create ext extrenally and all functions should use that directly
% - LQRcost_cov should take in the mean directly

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

rv_samples = 1000;
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
x_ol_stoch = zeros(nx, Tsim/Tfid + 1, rv_samples);
parfor i = 1:rv_samples
  [~, x_ode] = ode45(@(t, x) msd(t, x, 0, A_msd, B_msd), 0:Tfid:Tsim, x0_rv(:, i));
  x_ol_stoch(:, :, i) = x_ode';
end
data.x_ol_stoch = x_ol_stoch;

%% C.2 plot open loop dynamics
fig = figure;
hold on;
plot(data.t, squeeze(data.x_ol_stoch(1, :, :)), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.2]);
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
  x_cl_stoch = zeros(nx, Tsim/Tfid + 1, rv_samples);
  x0_rv_ext = [x0_rv; repmat(u0, 1, rv_samples); repmat(ref, 1, rv_samples)];
  data.cost_lqr = LQRCost_vec(x0_rv_ext, Uopt, Q_ext, S, M, Qbar, Rbar);
  parfor i = 1:rv_samples
    x_cl_stoch(:,:,i) = sim_closed_loop(Tsim, Tfid, Ts, x0_rv(:, i), u0, Uopt, A_msd, B_msd, @msd);
  end
  data.x_cl_stoch = x_cl_stoch;
  
  %% plot closed loop dynamics
  fig = figure;
  hold on;
  
  plot(data.t, squeeze(data.x_cl_stoch(1, :, :)), 'b', 'LineWidth', 1, 'HandleVisibility', 'off', 'Color', [0.5 0.5 0.5, 0.2]);
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
corr_st = calc_corr_st(cost_hf, cost_lf);
corr_an = calc_corr_an(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
% plot costs
title_str = ["$J_h$ and $J_l$", "$J_h$ (restriction)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hlr})$"];
plot_multifid_costs(perturb_range, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Perturbation");

%% E.1.1.2 averaging
Uopt_lf = downsample_avg(Uopt_hf, 10); % get HF in LF by avg every 10 values of HF
% vis_sols(Uopt_hf, Uopt_lf, data.lqrsol{1}.times, "HF obj", perturb_range, "Perturbation");
Uopt_lf = Uopt_lf(1:10:end, :); % pick every 10th elem
% calculate costs and correlation for HF/LF sols for each perturbation
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
corr_st = calc_corr_st(cost_hf, cost_lf);
corr_an = calc_corr_an(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
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
corr_st = calc_corr_st(cost_hf, cost_lf);
corr_an = calc_corr_an(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, Uopt_hf, Uopt_lf);
% plot costs
title_str = ["$J_h$ and $J_l$", "$J_l$"];
costs_str = ["$J_h(u_lh)$", "$J_l(u_l)$"];
plot_multifid_costs(perturb_range, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Perturbation");

%% *** E.2 correlation at different points in a numerical optimizer
%% E.2.1 num opt with HF objective fn
Uopt_hf = data.lqrsol{1}.Uopt;
u0_num = repelem(data.lqrsol{2}.Uopt, 10, 1); % warm start
% u0_num = zeros(size(Uopt_hf)); % zero start
% u0_num = randn(size(Uopt_hf)); % random start
global num_opt_iters num_opt_data num_opt_fvals;
num_opt_iters = 0;
num_opt_data = zeros(size(Uopt_hf,1), 100);
fun = @(u) LQRCostwithGrad(data.lqrsol{1}.x0_ext, u, data.lqrsol{1}.Q_ext, data.lqrsol{1}.S, data.lqrsol{1}.M, data.lqrsol{1}.Qbar, data.lqrsol{1}.Rbar);
options = optimoptions('fminunc', 'Display', 'iter', 'Algorithm', 'quasi-newton', 'SpecifyObjectiveGradient', true, 'OutputFcn', @outfun);
Uopt_num = fminunc(fun, u0_num, options);
uopt_diff = sqrt(sum((Uopt_hf - Uopt_num).^2));

U_hf = num_opt_data(:, 1:num_opt_iters);
U_num_hf = U_hf; % save for plotting later
% this is just to get the costs for the averaged HF solution in LF sim as the MLMC estimator would
%% E.2.1.1 correlation with restriction
U_lf = U_hf(1:10:end, :);
% U_lf = repelem(U_lf, 10, 1); % uncomment to visualize
% vis_sols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:num_opt_iters, "Iteration");
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
corr_st = calc_corr_st(cost_hf, cost_lf);
corr_an = calc_corr_an(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);

title_str = ["$J_h$ and $J_l$", "$J_h$ (restriction)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hlr})$"];
plot_multifid_costs(1:num_opt_iters, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Iteration");
% Get correlation for all iterations
corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
% plot correlation
plot_corr_2d(corr, "Correlation between costs in $J_h(u_h)$ and $J_l(u_{hlr})$", costs_str, "num_iters_res");

%% E.2.1.2 correlation with averaging
U_lf = downsample_avg(U_hf, 10); % get LF by avg every 10 values of HF
% vis_sols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:num_opt_iters, "Iteration");
U_lf = U_lf(1:10:end, :);
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
corr_st = calc_corr_st(cost_hf, cost_lf);
corr_an = calc_corr_an(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);

title_str = ["$J_h$ and $J_l$", "$J_h$ (averaging)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla})$"];
plot_multifid_costs(1:num_opt_iters, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Iteration");
corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
% plot correlation
plot_corr_2d(corr, "Correlation between costs in $J_h(u_h)$ and $J_l(u_{hla})$", costs_str, "num_iters_avg");

%% F num opt with MLMC estimator
%% F.1 normal MLMC Estimator
cv_samples = 500;
num_opt_iters = 0;
num_opt_data = zeros(size(Uopt_hf,1), 100);
global alpha_mc;
alpha_mc = zeros(1, 100);
fun = @(u) calc_mlmc_est(data.lqrsol{1}, data.lqrsol{2}, u, cv_samples, cv_samples, u0, ref, x0_rv, x0_mean, x0_cov);
options = optimoptions('fminunc', 'Display', 'iter', 'OutputFcn', @outfun); % finite diff
Uopt_num = fminunc(fun, u0_num, options);
U_hf = num_opt_data(:, 1:num_opt_iters);
U_num_mlmc = U_hf; % save for plotting later
U_lf = downsample_avg(U_hf, 10);
% vis_sols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:num_opt_iters, "Iteration");

U_lf = U_lf(1:10:end, :);
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
corr_st = calc_corr_st(cost_hf, cost_lf);
corr_an = calc_corr_an(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);

title_str = ["$J_h$ and $J_l$", "$S_{" + cv_samples + "}^{CV}$"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla})$"];
plot_multifid_costs(1:num_opt_iters, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Iteration");
corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
% plot correlation
plot_corr_2d(corr, "$S_{" + cv_samples + "}^{CV}$ opt Correlation between costs in $J_h(u_h)$ and $J_l(u_{hla})$", costs_str, "num_iters_mlmc");

figure;
plot(alpha_mc(1:num_opt_iters), 'LineWidth', 2);
title("MLMC estimator alpha values");
xlabel("Iteration");
ylabel("Alpha");
grid on;

%% F.2 MLMC Estimator with fixed LF solution
[~, it_max_cor] = max(sum(corr, 1));
u_lf = U_lf(:, it_max_cor);
num_opt_iters = 0;
num_opt_data = zeros(size(Uopt_hf,1), 100);
alpha_mc = zeros(1, 100);
fun = @(u) calc_mlmc_est_max_corr(data.lqrsol{1}, data.lqrsol{2}, u, u_lf, cv_samples, cv_samples, u0, ref, x0_rv, x0_mean, x0_cov);
options = optimoptions('fminunc', 'Display', 'iter', 'OutputFcn', @outfun); % finite dif
Uopt_num = fminunc(fun, u0_num, options);
U_hf = num_opt_data(:, 1:num_opt_iters);
U_num_mlmc_fix = U_hf; % save for plotting later
U_lf = repelem(u_lf, 10, num_opt_iters);
% vis_sols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:num_opt_iters, "Iteration");

U_lf = U_lf(1:10:end, :);
[cost_hf, cost_lf] = calc_costs_multifid(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
corr_st = calc_corr_st(cost_hf, cost_lf);
corr_an = calc_corr_an(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);

title_str = ["$J_h$ and $J_l$ at max correlation", "$S_{" + cv_samples + "}^{CV}$"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla}^{max})$"];
plot_multifid_costs(1:num_opt_iters, mean(cost_hf,2), mean(cost_lf,2), corr_st, corr_an, title_str, costs_str, "Iteration");
corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, data.lqrsol{1}, data.lqrsol{2}, U_hf, U_lf);
% plot correlation
plot_corr_2d(corr, "$S_{" + cv_samples + "}^{CV}$ opt Correlation between costs in $J_h(u_h)$ and $J_l(u_{hla}^{max})$", costs_str, "num_iters_mlmc_max_corr");

figure;
plot(alpha_mc(1:num_opt_iters), 'LineWidth', 2);
title("MLMC estimator at max corelation alpha values");
xlabel("Iteration");
ylabel("Alpha");
grid on;

%% plot convergence distance comparison
plot_convergence_dist(data.lqrsol{1}.Uopt, U_num_hf, U_num_mlmc, U_num_mlmc_fix, ["HF", "CV", "CV with max corr"], "Convergence distance comparison");

%% G Convergence and variance with various optimizers and sample sizes
num_rv_samples = [100];
num_estimator_samples = 100;
u0_num = repelem(data.lqrsol{2}.Uopt, 10, 1); % warm start
max_iters = 30;
options = optimoptions('fminunc', 'Display', 'iter', 'OutputFcn', @outfun, "MaxIterations", max_iters);
data.hf_var = zeros(length(num_rv_samples), max_iters);
data.mlmc_var = zeros(length(num_rv_samples), max_iters);
data.mlmc_fix_var = zeros(length(num_rv_samples), max_iters);
data.hf_u = zeros(length(u0_num), max_iters, num_estimator_samples, length(num_rv_samples));
data.mlmc_u = zeros(length(u0_num), max_iters, num_estimator_samples, length(num_rv_samples));
data.mlmc_fix_u = zeros(length(u0_num), max_iters, num_estimator_samples, length(num_rv_samples));
for num_samples=num_rv_samples
  fprintf("\n*** RV Samples: %d ***\n", num_samples);
  wait_bar = waitbar(0, "Running estimators with " + num_samples + " samples");
  temp_cost.hf = zeros(num_estimator_samples, max_iters);
  temp_cost.mlmc = zeros(num_estimator_samples, max_iters);
  temp_cost.mlmc_fix = zeros(num_estimator_samples, max_iters);
  for i=1:num_estimator_samples
    waitbar(i/num_estimator_samples, wait_bar);
    % get the RV samples
    x0_rv = mvnrnd(x0_mean, x0_cov, num_samples)';
    x0_rv_ext = [x0_rv; repmat(u0, 1, num_samples); repmat(ref, 1, num_samples)];
    
    % for just the HF cost fn
    num_opt_iters = 0;
    num_opt_data = zeros(size(Uopt_hf,1), max_iters);
    num_opt_fvals = zeros(1, max_iters);
    x0_rv_ext_mean = mean(x0_rv_ext, 2);
    fun = @(u) mean(LQRCost_vec(x0_rv_ext, u, data.lqrsol{1}.Q_ext, data.lqrsol{1}.S, data.lqrsol{1}.M, data.lqrsol{1}.Qbar, data.lqrsol{1}.Rbar));
    Uopt_num = fminunc(fun, u0_num, options);
    temp_cost.hf(i, :) = num_opt_fvals(1:max_iters);
    data.hf_u(:, :, i, num_rv_samples == num_samples) = num_opt_data(:, 1:max_iters);
    
    % for MLMC estimator
    num_opt_iters = 0;
    num_opt_data = zeros(size(Uopt_hf,1), max_iters);
    num_opt_fvals = zeros(1, max_iters);
    fun = @(u) calc_mlmc_est(data.lqrsol{1}, data.lqrsol{2}, u, num_samples, num_samples, u0, ref, x0_rv, x0_mean, x0_cov);
    Uopt_num = fminunc(fun, u0_num, options);
    temp_cost.mlmc(i, :) = num_opt_fvals(1:max_iters);
    data.mlmc_u(:, :, i, num_rv_samples == num_samples) = num_opt_data(:, 1:max_iters);
    
    % for MLMC estimator with fixed LF solution
    u_lf = num_opt_data(:, it_max_cor); % take from prev
    u_lf = u_lf(1:10:end, :);
    num_opt_iters = 0;
    num_opt_data = zeros(size(Uopt_hf,1), max_iters);
    num_opt_fvals = zeros(1, max_iters);
    fun = @(u) calc_mlmc_est_max_corr(data.lqrsol{1}, data.lqrsol{2}, u, u_lf, num_samples, num_samples, u0, ref, x0_rv, x0_mean, x0_cov);
    Uopt_num = fminunc(fun, u0_num, options);
    temp_cost.mlmc_fix(i, :) = num_opt_fvals(1:max_iters);
    data.mlmc_fix_u(:, :, i, num_rv_samples == num_samples) = num_opt_data(:, 1:max_iters);
  end
  data.hf_var(num_rv_samples == num_samples, :) = var(temp_cost.hf, 0, 1);
  data.mlmc_var(num_rv_samples == num_samples, :) = var(temp_cost.mlmc, 0, 1);
  data.mlmc_fix_var(num_rv_samples == num_samples, :) = var(temp_cost.mlmc_fix, 0, 1);
  close(wait_bar);
end

%% calculate analytical variance
data.hf_var_an = zeros(length(num_rv_samples), max_iters);
data.mlmc_var_an = zeros(length(num_rv_samples), max_iters);
data.mlmc_fix_var_an = zeros(length(num_rv_samples), max_iters);
x0_rv_mean_ext = [x0_mean; u0; ref];
x0_rv_cov_ext = blkdiag(x0_cov, zeros(3,3));
for i=1:length(num_rv_samples)
  for j=1:max_iters
    % average out u vlues for calculating variance
    hf_u_avg = mean(data.hf_u(:, j, :, i), 3);
    data.hf_var_an(i, j) = mc_var(data.lqrsol{1}, x0_rv_mean_ext, x0_rv_cov_ext, hf_u_avg, num_rv_samples(i));
    
    mlmc_u_avg = mean(data.mlmc_u(:, j, :, i), 3);
    data.mlmc_var_an(i, j) = cv_var(data.lqrsol{1}, data.lqrsol{2}, x0_rv_ext, x0_rv_cov_ext, mlmc_u_avg, num_rv_samples(i));
    
    u_lf_avg = mean(data.mlmc_u(:, it_max_cor, :, i), 3);
    u_lf_avg = u_lf_avg(1:10:end, :);
    mlmc_fix_u_avg = mean(data.mlmc_fix_u(:, j, :, i), 3);
    data.mlmc_fix_var_an(i, j) = cv_max_var(data.lqrsol{1}, data.lqrsol{2}, x0_rv_ext, x0_rv_cov_ext, mlmc_fix_u_avg, u_lf_avg, num_rv_samples(i));
  end
end

%% plot variance
for i=1:length(num_rv_samples)
  fig = figure;
  semilogy(1:max_iters, data.hf_var(i, :), 'b', 'LineWidth', 2, 'DisplayName', 'HF');
  hold on;
  semilogy(1:max_iters, data.mlmc_var(i, :), 'r', 'LineWidth', 2, 'DisplayName', 'CV');
  semilogy(1:max_iters, data.mlmc_fix_var(i, :), 'g', 'LineWidth', 2, 'DisplayName', 'CV with max corr');
  semilogy(1:max_iters, data.hf_var_an(i, :), '--b', 'LineWidth', 2, 'DisplayName', 'HF (Analytical)');
  semilogy(1:max_iters, data.mlmc_var_an(i, :), '--r', 'LineWidth', 2, 'DisplayName', 'CV (Analytical)');
  semilogy(1:max_iters, data.mlmc_fix_var_an(i, :), '--g', 'LineWidth', 2, 'DisplayName', 'CV with max corr (Analytical)');
  title("Cost variance vs iterations with "+num_rv_samples(i)+" samples");
  xlabel("Iteration");
  ylabel("Variance");
  legend show;
  grid on;
end

%% plot convergence in solution
for i=1:length(num_rv_samples)
  plot_convergence_dist(data.lqrsol{1}.Uopt, squeeze(data.hf_u(:, :, :, i)), squeeze(data.mlmc_u(:, :, :, i)), squeeze(data.mlmc_fix_u(:, :, :, i)), ...
    ["HF", "CV", "CV with max corr"], "Convergence distance comparison with " + num_rv_samples(i) + "samples");
end

%% plot convergence in cost

%% save all data
save("artefacts/data.mat");

function xdot = msd(~, x, u, A, B)
xdot = A * x + B * u;
end

function cost = LQRCost(x0, u, Q, S, M, Qbar, Rbar)
cost = u'*(S'*Qbar*S + Rbar)*u + 2*x0'*M'*Qbar*S*u + x0'*(M'*Qbar*M + Q)*x0;
end

function cost = LQRCost_vec(x0, u, Q, S, M, Qbar, Rbar) % expects x0 to be multiple samples
cost = u'*(S'*Qbar*S + Rbar)*u + 2*x0'*M'*Qbar*S*u + diag(x0'*(M'*Qbar*M + Q)*x0);
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

function stop = outfun(x, optimValues, state)
% this is for capturing intermediate values of the optimization
stop = false;
global num_opt_iters num_opt_data num_opt_fvals;
if isequal(state, 'iter')
  num_opt_iters = num_opt_iters + 1;
  num_opt_data(:, num_opt_iters) = x;
  num_opt_fvals(num_opt_iters) = optimValues.fval;
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

function corr_st = calc_corr_st(cost_hf, cost_lf)
perturbs = size(cost_hf, 1);
corr_st = zeros(perturbs, 1);
for i = 1:perturbs
  corr_mat = corrcoef(cost_hf(i, :), cost_lf(i, :));
  corr_st(i) = corr_mat(1,2); % we only need the cross correlation, diagnonal will be 1
end
end

function corr_an = calc_corr_an(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, U_hf, U_lf)
perturbs = size(U_hf, 2);
corr_an = zeros(perturbs, 1);
x0_rv_ext = [x0_rv; repmat(u0, 1, size(x0_rv, 2)); repmat(ref, 1, size(x0_rv, 2))];
x_mean = mean(x0_rv_ext, 2);
x_cov = cov(x0_rv_ext');
for i = 1:perturbs
  cov_an = LQRcost_cov(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, U_hf(:, i), U_lf(:, i)); % TODO make indices consistent
  [~, var_J1] = LQRcost_stats(x_mean, x_cov, U_hf(:, i), lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
  [~, var_J2] = LQRcost_stats(x_mean, x_cov, U_lf(:, i), lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
  corr_an(i) = cov_an / sqrt(var_J1 * var_J2);
end
end


function corr = calc_corr_multifid_2d_iters(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, Uopt_hf, Uopt_lf)
iters = size(Uopt_hf, 2);
rv_samples = size(x0_rv, 2);
cost_hf_corr = zeros(iters, rv_samples);
cost_lf_corr = zeros(iters, rv_samples);
x0_rv_ext = [x0_rv; repmat(u0, 1, rv_samples); repmat(ref, 1, rv_samples)];
for i = 1:iters
  cost_hf_corr(i, :) = LQRCost_vec(x0_rv_ext, Uopt_hf(:, i), lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
  cost_lf_corr(i, :) = LQRCost_vec(x0_rv_ext, Uopt_lf(:, i), lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
end

corr = zeros(iters, iters);
for i = 1:iters % HF iterations
  for j = 1:iters % LF iterations
    corr_mat = corrcoef(cost_hf_corr(i, :), cost_lf_corr(j, :));
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
x0_rv_ext = [x0_rv; repmat(u0, 1, rv_samples); repmat(ref, 1, rv_samples)];
for i = 1:perturbs
  cost_hf(i,:) = LQRCost_vec(x0_rv_ext, Uopt_hf(:, i), lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
  cost_lf(i,:) = LQRCost_vec(x0_rv_ext, Uopt_lf(:, i), lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
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

function cost = calc_mlmc_est(lqrsol_hf, lqrsol_lf, u_hf, n_hf, n_lf, u0, ref, x0_rv, x0_rv_mean, x0_rv_cov)
% MC estimator for hf
x0_rv_ext = [x0_rv(:, 1:n_hf); repmat(u0, 1, n_hf); repmat(ref, 1, n_hf)];
cost_hf = mean(LQRCost_vec(x0_rv_ext, u_hf, lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar));

% reuse same samples for LF
u_lf = downsample_avg(u_hf, 10);
u_lf = u_lf(1:10:end);
cost_lf = mean(LQRCost_vec(x0_rv_ext, u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar));

% use all samples for expectation
x0_rv_mean_ext = [x0_rv_mean; u0; ref];
x0_rv_cov_ext = [x0_rv_cov, zeros(size(x0_rv_cov, 1), size(u0, 1) + size(ref, 1));
  zeros(size(u0, 1) + size(ref, 1), size(x0_rv_cov, 1) + size(u0, 1) + size(ref, 1))];
[exp_l, var_l] = LQRcost_stats(x0_rv_mean_ext, x0_rv_cov_ext, u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);

% calculate optimal alpha
cov_hl = LQRcost_cov(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, u_hf, u_lf);
alpha = -cov_hl / var_l;
global alpha_mc num_opt_iters;
alpha_mc(num_opt_iters+1) = alpha;
cost = cost_hf + alpha * (cost_lf - exp_l);

end

function cost = calc_mlmc_est_max_corr(lqrsol_hf, lqrsol_lf, u_hf, u_lf, n_hf, n_lf, u0, ref, x0_rv, x0_rv_mean, x0_rv_cov)
% MC estimator for hf
x0_rv_ext = [x0_rv(:, 1:n_hf); repmat(u0, 1, n_hf); repmat(ref, 1, n_hf)];
cost_hf = LQRCost_vec(x0_rv_ext, u_hf, lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
cost_hf = mean(cost_hf);

% reuse same samples for LF
% u_lf = downsample_avg(u_hf, 10);
% u_lf = u_lf(1:10:end);
cost_lf = LQRCost_vec(x0_rv_ext, u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);
cost_lf = mean(cost_lf);

% use all samples for expectation
x0_rv_mean_ext = [x0_rv_mean; u0; ref];
x0_rv_cov_ext = [x0_rv_cov, zeros(size(x0_rv_cov, 1), size(u0, 1) + size(ref, 1));
  zeros(size(u0, 1) + size(ref, 1), size(x0_rv_cov, 1) + size(u0, 1) + size(ref, 1))];
[exp_l, var_l] = LQRcost_stats(x0_rv_mean_ext, x0_rv_cov_ext, u_lf, lqrsol_lf.Q_ext, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar, lqrsol_lf.Rbar);


% calculate optimal alpha
cov_hl = LQRcost_cov(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, u_hf, u_lf);
alpha = -cov_hl / var_l;
global alpha_mc num_opt_iters;
alpha_mc(num_opt_iters+1) = alpha;
cost = cost_hf + alpha * (cost_lf - exp_l);

end

function plot_convergence_dist(an_sol, U_num_hf, U_num_mlmc, U_num_mlmc_fix, labels, title_str)
dist_hf = squeeze(vecnorm(an_sol - U_num_hf, 2, 1));
dist_mlmc = squeeze(vecnorm(an_sol - U_num_mlmc, 2, 1));
dist_mlmc_fix = squeeze(vecnorm(an_sol - U_num_mlmc_fix, 2, 1));
if size(dist_hf, 1) == 1
  dist_hf_mean = dist_hf;
  dist_mlmc_mean = dist_mlmc;
  dist_mlmc_fix_mean = dist_mlmc_fix;
else
  dist_hf_mean = mean(dist_hf, 2);
  dist_mlmc_mean = mean(dist_mlmc, 2);
  dist_mlmc_fix_mean = mean(dist_mlmc_fix, 2);
end
figure;
% handle alpha to be min 0.1 and max 1 based on how many lines there are
semilogy(dist_hf, 'LineWidth', 1, 'Color', [0 0 1, max(0.1, 1/size(dist_hf, 1))]);
hold on;
% plot means with thicker lines
semilogy(dist_mlmc, 'LineWidth', 1, 'Color', [1 0 0, max(0.1, 1/size(dist_mlmc, 1))]);
semilogy(dist_mlmc_fix, 'LineWidth', 1, 'Color', [0 1 0, max(0.1, 1/size(dist_mlmc_fix, 1))]);
l{1} = semilogy(dist_hf_mean, 'b', "Marker", "x", 'LineWidth', 2);
l{2} = semilogy(dist_mlmc_mean, 'r', "Marker", "x", 'LineWidth', 2);
l{3} = semilogy(dist_mlmc_fix_mean, 'g', "Marker", "x", 'LineWidth', 2);
title(title_str);
xlabel("Iteration");
ylabel("Distance to analytical solution");
legend([l{:}], labels);
grid on;
end

function var = mc_var(lqrsol, x0_rv_ext_mean, x0_rv_ext_cov, u, num_samples)
[~, var] = LQRcost_stats(x0_rv_ext_mean, x0_rv_ext_cov, u, lqrsol.Q_ext, lqrsol.S, lqrsol.M, lqrsol.Qbar, lqrsol.Rbar);
var = var / num_samples;
end

function var = cv_var(lqrsol_hf, lqrsol_lf, x0_rv_ext_mean, x0_rv_ext_cov, u, n_hf)
wtf = mean(x0_rv_ext_mean, 2);
[~, var_hf] = LQRcost_stats(wtf, x0_rv_ext_cov, u, lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
u_lf = downsample_avg(u, 10);
u_lf = u_lf(1:10:end);
x0_rv_mean = x0_rv_ext_mean(1:2, :);
u0 = x0_rv_ext_mean(3, 1);
ref = x0_rv_ext_mean(4:5, 1);
corr = calc_corr_an(x0_rv_mean, u0, ref, lqrsol_hf, lqrsol_lf, u, u_lf);
var = var_hf/n_hf * (1 - corr^2);
end

function var = cv_max_var(lqrsol_hf, lqrsol_lf, x0_rv_ext_mean, x0_rv_ext_cov, u, u_lf, n_hf)
wtf = mean(x0_rv_ext_mean, 2);
[~, var_hf] = LQRcost_stats(wtf, x0_rv_ext_cov, u, lqrsol_hf.Q_ext, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar, lqrsol_hf.Rbar);
x0_rv_mean = x0_rv_ext_mean(1:2, :);
u0 = x0_rv_ext_mean(3, 1);
ref = x0_rv_ext_mean(4:5, 1);
corr = calc_corr_an(x0_rv_mean, u0, ref, lqrsol_hf, lqrsol_lf, u, u_lf);
var = var_hf/n_hf * (1 - corr^2);
end
