% TODOS:
% - we are assuming u0 to be 0, might not be true in general so might need to add u0 to cumsum
% - make indices consistent when storing U in row vs column (make everything row based)
% - calculate cost of CV estimators to determine equal cost split

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
u0 = 0; % initial control effort !=0 to prevent repmat warnings
x0_ext = [x0; u0; ref]; % extended state
x0_ext_mean = [x0_mean; u0; ref];
x0_ext_cov = blkdiag(x0_cov, zeros(length(u0), length(u0)), zeros(length(ref), length(ref)));

Tsim = 6;
Tfid = 0.01; % simulation fidelity
% Ts = 0.05; % sampling time, MUST BE MULTIPLE OF Tfid
% Ts = 0.5; % sampling time, MUST BE MULTIPLE OF Tfid
TsList = [0.05, 0.5];

Qe = eye(nx) * 2; % state cost
Re = eye(nu) * 10; % input cost
% terminal state cost, needs to be 0 for the extended state since
% we dont want to penalize the state being away from the origin
Pe = eye(nx) * 0;

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
x0_rv_ext = [x0_rv; repmat(u0, 1, rv_samples); repmat(ref, 1, rv_samples)];
x0_rv_ext_mean = [x0_rv_mean; u0; ref];
x0_rv_ext_cov = blkdiag(x0_rv_cov, zeros(length(u0), length(u0)), zeros(length(ref), length(ref)));

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
  [A_ext, B_ext, Q, R, P] = extendState(c2d(msd_sys, Ts), Qe, Re, Pe);
  [Kopt, S, M, Qbar, Rbar] = solveLQR(N, A_ext, B_ext, Q, R, P);
  Uopt = Kopt * x0_ext;
  cost_det = St.LQRCost(x0_ext, struct('Q', Q, 'S', S, 'M', M, 'Qbar', Qbar, 'Rbar', Rbar), Uopt);
  fprintf("Deterministic LQR cost: %f\n", cost_det);
  
  %% simulate closed loop dynamics
  data.x_cl = simCL(Tsim, Tfid, Ts, x0, u0, Uopt, A_msd, B_msd, @msd);
  
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
  Uopt = Kopt * x0_rv_ext_mean;
  data.lqrsol{idx} = struct('x0_ext', x0_rv_ext_mean, 'Uopt', Uopt, 'Q', Q, 'S', S, 'M', M, 'Qbar', Qbar, 'Rbar', Rbar, 'times', times);
  x_cl_stoch = zeros(nx, Tsim/Tfid + 1, rv_samples);
  data.cost_lqr = St.LQRCost(x0_rv_ext, data.lqrsol{idx}, Uopt);
  parfor i = 1:rv_samples
    x_cl_stoch(:,:,i) = simCL(Tsim, Tfid, Ts, x0_rv(:, i), u0, Uopt, A_msd, B_msd, @msd);
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
  cost_lqr_exp = St.LQRExp(x0_rv_ext_mean, x0_rv_ext_cov, data.lqrsol{idx}, Uopt);
  fprintf("Expectaion of Stochastic LQR cost\n- Analytical: %f\n- Experimental: %f\n", cost_lqr_exp, mean(data.cost_lqr));
  cost_lqr_var = St.LQRVar(x0_rv_ext_mean, x0_rv_ext_cov, data.lqrsol{idx}, Uopt);
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

% *** E. correlation check for CV
% E.1 correlation with perturbation in direction of maximum ascent

%% E.1.1 perturbation in HF objective fn
perturbation = perturbMaxGrad(data.lqrsol{1}, perturb_dir_samples, perturb_dir_mag, perturb_range);
U_h = data.lqrsol{1}.Uopt + perturbation; % perturb hf sol

%% E.1.1.1 restriction
U_l = U_h(1:10:end, :); % pick every 10th elem
% U_l_vis = repelem(U_l, 10, 1); % uncomment to visualize
% visSols(U_h, U_l_vis , data.lqrsol{1}.times, "HF obj", perturb_range, "Perturbation");

% calculate costs and correlation for HF/LF sols for each perturbation
cost_h = St.LQRCostMulti(x0_rv_ext, data.lqrsol{1}, U_h);
cost_l = St.LQRCostMulti(x0_rv_ext, data.lqrsol{2}, U_l);
corr_st = St.CorrMulti(cost_h, cost_l);
corr_an = St.LQRCorrMulti(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);
% plot costs
title_str = ["$J_h$ and $J_l$", "$J_h$ (restriction)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hlr})$"];
plotMultifidCost(perturb_range, mean(cost_h,2), mean(cost_l,2), corr_st, corr_an, title_str, costs_str, "Perturbation");

%% E.1.1.2 averaging
U_l = Est.DownsampleAvg(U_h, 10); % get HF in LF by avg every 10 values of HF
% Uopt_lf = Est.DownsampleAvg(Uopt_hf, 10, true); % uncomment to visualize
% visSols(Uopt_hf, Uopt_lf, data.lqrsol{1}.times, "HF obj", perturb_range, "Perturbation");
% calculate costs and correlation for HF/LF sols for each perturbation
cost_h = St.LQRCostMulti(x0_rv_ext, data.lqrsol{1}, U_h);
cost_l = St.LQRCostMulti(x0_rv_ext, data.lqrsol{2}, U_l);
corr_st = St.CorrMulti(cost_h, cost_l);
corr_an = St.LQRCorrMulti(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);
% plot costs
title_str = ["$J_h$ and $J_l$", "$J_h (averaging)$"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla})$"];
plotMultifidCost(perturb_range, mean(cost_h,2), mean(cost_l,2), corr_st, corr_an, title_str, costs_str, "Perturbation");


%% E.1.2 perturbation in LF objective fn
perturbation = perturbMaxGrad(data.lqrsol{2}, perturb_dir_samples, perturb_dir_mag, perturb_range);
U_l = data.lqrsol{2}.Uopt + perturbation;
% visSols(Uopt_hf, Uopt_lf, data.lqrsol{2}.times, "LF obj (restriction)", perturb_range, "Perturbation");
U_h = repelem(U_l, 10, 1); % repeat the LF sol to match HF sol
% calculate costs and correlation for both high and low fidelity, for each perturbation
cost_h = St.LQRCostMulti(x0_rv_ext, data.lqrsol{1}, U_h);
cost_l = St.LQRCostMulti(x0_rv_ext, data.lqrsol{2}, U_l);
corr_st = St.CorrMulti(cost_h, cost_l);
corr_an = St.LQRCorrMulti(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);
% plot costs
title_str = ["$J_h$ and $J_l$", "$J_l$"];
costs_str = ["$J_h(u_lh)$", "$J_l(u_l)$"];
plotMultifidCost(perturb_range, mean(cost_h,2), mean(cost_l,2), corr_st, corr_an, title_str, costs_str, "Perturbation");

%% *** E.2 correlation at different points in a numerical optimizer
%% E.2.1 num opt with HF objective fn
u0_num = repelem(data.lqrsol{2}.Uopt, 10, 1); % warm start
fun = @(u) St.LQRCostwGrad(data.lqrsol{1}.x0_ext, data.lqrsol{1}, u);
U_num = fminuncWHistory(fun, u0_num);
Uopt_h_num = U_num(:, end);
iters = size(U_num, 2);

uopt_diff = mean((data.lqrsol{1}.Uopt - Uopt_h_num).^2);

% this is just to get the costs for the averaged HF solution in LF sim as the CV estimator would
%% E.2.1.1 correlation with restriction
U_h = U_num;
U_l = U_h(1:10:end, :);
% U_lf = repelem(U_lf, 10, 1); % uncomment to visualize
% visSols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:iters, "Iteration");
cost_h = St.LQRCostMulti(x0_rv_ext, data.lqrsol{1}, U_h);
cost_l = St.LQRCostMulti(x0_rv_ext, data.lqrsol{2}, U_l);
corr_st = St.CorrMulti(cost_h, cost_l);
corr_an = St.LQRCorrMulti(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);
corr_2d = St.LQRCorrMulti2D(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);

% plot correlation
title_str = ["$J_h$ and $J_l$", "$J_h$ (restriction)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hlr})$"];
plotMultifidCost(1:iters, mean(cost_h,2), mean(cost_l,2), corr_st, corr_an, title_str, costs_str, "Iteration");
plotCorr2d(corr_2d, "Correlation between costs in $J_h(u_h)$ and $J_l(u_{hlr})$", costs_str, "num_iters_res");

%% E.2.1.2 correlation with averaging
U_l = Est.DownsampleAvg(U_h, 10); % get LF by avg every 10 values of HF
% U_lf = Est.DownsampleAvg(U_hf, 10, true); % uncomment to visualize
% visSols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:iters, "Iteration");
cost_h = St.LQRCostMulti(x0_rv_ext, data.lqrsol{1}, U_h);
cost_l = St.LQRCostMulti(x0_rv_ext, data.lqrsol{2}, U_l);
corr_st = St.CorrMulti(cost_h, cost_l);
corr_an = St.LQRCorrMulti(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);
corr_2d = St.LQRCorrMulti2D(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);

% plot correlation
title_str = ["$J_h$ and $J_l$", "$J_h$ (averaging)"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla})$"];
plotMultifidCost(1:iters, mean(cost_h,2), mean(cost_l,2), corr_st, corr_an, title_str, costs_str, "Iteration");
plotCorr2d(corr_2d, "Correlation between costs in $J_h(u_h)$ and $J_l(u_{hla})$", costs_str, "num_iters_avg");

%% F num opt with CV estimator
cv = Cv(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2});
%% F.1 normal CV Estimator
cv_samples = rv_samples/2;
[~, U_h, U_l] = cv.opt(u0_num, -1, -1, x0_rv_ext, cv_samples, false);
U_num_cv = U_h; % save for plotting later

% U_lf = repelem(U_hla, 10, cv.idx); % uncomment to visualize
% visSols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:iters, "Iteration");

cost_h = St.LQRCostMulti(x0_rv_ext, data.lqrsol{1}, U_h);
cost_l = St.LQRCostMulti(x0_rv_ext, data.lqrsol{2}, U_l);
corr_st = St.CorrMulti(cost_h, cost_l);
corr_an = St.LQRCorrMulti(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);
corr_2d = St.LQRCorrMulti2D(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);

% plot correlation
title_str = ["$J_h$ and $J_l$", "$S_{" + cv_samples + "}^{CV}$"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla})$"];
plotMultifidCost(1:cv.idx, mean(cost_h,2), mean(cost_l,2), corr_st, corr_an, title_str, costs_str, "Iteration");
plotCorr2d(corr_2d, "$S_{" + cv_samples + "}^{CV}$ opt Correlation between costs in $J_h(u_h)$ and $J_l(u_{hla})$", costs_str, "num_iters_cv");

%% F.2 CV Estimator with fixed LF solution
[~, U_h, U_l] = cv.opt(u0_num, -1, -1, x0_rv_ext, cv_samples, true);
U_num_cv_fix = U_h; % save for plotting later

% U_lf = repelem(U_hla, 10, cv.idx);
% visSols(U_hf, U_lf, data.lqrsol{1}.times, "Solutions along optimizer path", 1:iters, "Iteration");

cost_h = St.LQRCostMulti(x0_rv_ext, data.lqrsol{1}, U_h);
cost_l = St.LQRCostMulti(x0_rv_ext, data.lqrsol{2}, U_l);
corr_st = St.CorrMulti(cost_h, cost_l);
corr_an = St.LQRCorrMulti(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);
corr_2d = St.LQRCorrMulti2D(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2}, U_h, U_l);

% plot correlation
title_str = ["$J_h$ and $J_l$ at max correlation", "$S_{" + cv_samples + "}^{CV}$"];
costs_str = ["$J_h(u_h)$", "$J_l(u_{hla}^{max})$"];
plotMultifidCost(1:cv.idx, mean(cost_h,2), mean(cost_l,2), corr_st, corr_an, title_str, costs_str, "Iteration");
plotCorr2d(corr_2d, "$S_{" + cv_samples + "}^{CV}$ opt Correlation between costs in $J_h(u_h)$ and $J_l(u_{hla}^{max})$", costs_str, "num_iters_cv_max_corr");

%% plot convergence distance comparison
plotConvergence(data.lqrsol{1}.Uopt, U_num, U_num_cv, U_num_cv_fix, ["HF", "CV", "CV with max corr"], "Convergence distance comparison");

%% G verify ACV estimator
cv = Cv(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2});
acv = Acv(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2});

%% check variance of ACV across different ratios of HF/LF cost
acv_ratios = 0.1:0.1:2;
num_rv_samples = 2000; % large enough to give samples to acv
num_estimator_samples = 500;
x0_rv_ext = zeros(length(x0_ext), num_rv_samples, num_estimator_samples);
for i=1:num_estimator_samples
  x0_rv = mvnrnd(x0_mean, x0_cov, num_rv_samples)';
  x0_rv_ext(:, :, i) = [x0_rv; repmat(u0, 1, num_rv_samples); repmat(ref, 1, num_rv_samples)];
end

n_mc = 500;
lf_hf_cost_ratio = 0.2;
var_data_acv = zeros(length(acv_ratios), num_estimator_samples);
for i=1:length(acv_ratios)
  fprintf("\n*** ACV Ratio: %f ***\n", acv_ratios(i));
  n_acv = round(n_mc / (1 + lf_hf_cost_ratio + acv_ratios(i)*lf_hf_cost_ratio));
  m_acv = round(acv_ratios(i) * n_acv);
  for j=1:num_estimator_samples
    var_data_acv(i, j) = acv.est(x0_rv_ext(:, :, j), n_acv, m_acv, Uopt_h_num, false, '-1', 'stat');
  end
end

var_acv = var(var_data_acv, 0, 2);

%%
figure;
plot(acv_ratios, var_acv, 'b', 'LineWidth', 2, 'DisplayName', 'ACV');
xlabel("ACV Ratio (m:n)");
ylabel("Variance");
title("ACV estimator Variance across m:n ratios (eq cost, $\alpha=-1$)", "Interpreter", "latex");
legend show;
grid on;

%% plot variance of CV vs ACV with n_acv fixed and increasing m_acv
m_acvs = 100:100:2000;
num_rv_samples = max(m_acvs)*2;
num_estimator_samples = 500;
x0_rv_ext = zeros(length(x0_ext), num_rv_samples, num_estimator_samples);
for i=1:num_estimator_samples
  x0_rv = mvnrnd(x0_mean, x0_cov, num_rv_samples)';
  x0_rv_ext(:, :, i) = [x0_rv; repmat(u0, 1, num_rv_samples); repmat(ref, 1, num_rv_samples)];
end

n_mc = 500;
lf_hf_cost_ratio = 0.2;
n_cv = round(n_mc / (1 + lf_hf_cost_ratio));
n_acv = n_cv;
var_data_cv = zeros(length(m_acvs), num_estimator_samples);
var_data_acv = zeros(length(m_acvs), num_estimator_samples);
for i=1:length(m_acvs)
  fprintf("\n*** ACV m: %d ***\n", m_acvs(i));
  for j=1:num_estimator_samples
    var_data_cv(i, j) = cv.est(x0_rv_ext(:, :, j), n_cv, Uopt_h_num, false);
    var_data_acv(i, j) = acv.est(x0_rv_ext(:, :, j), n_acv, m_acvs(i), Uopt_h_num, false, '-1', 'stat');
  end
end

var_acv = var(var_data_acv, 0, 2);
var_cv = var(var_data_cv, 0, 2);
var_diff = (var_acv - var_cv).^2;

%%
figure;
semilogy(m_acvs, var_diff, 'r', 'LineWidth', 2, 'DisplayName', '(var[ACV] - var[CV])^2');
xlabel("m_{acv}");
ylabel("Variance");
title("ACV-CV Variance difference with different $m_{acv} (\alpha=-1)$", "Interpreter", "latex");
legend show;
grid on;

%% plot variance of CV/ACV vs analyitical
num_estimator_samples = 100:100:2000;
lf_hf_cost_ratio = 0.2;
acv_ratio  = 1;
n_mc = 500;
n_cv = round(n_mc / (1 + lf_hf_cost_ratio));
n_acv = round(n_mc / (1 + lf_hf_cost_ratio + acv_ratio*lf_hf_cost_ratio));
m_acv = round(acv_ratio * n_acv);
num_rv_samples = round(max(n_acv+m_acv, n_mc));
x0_rv_ext = zeros(length(x0_ext), num_rv_samples, max(num_estimator_samples));
for i=1:max(num_estimator_samples)
  x0_rv = mvnrnd(x0_mean, x0_cov, num_rv_samples)';
  x0_rv_ext(:, :, i) = [x0_rv; repmat(u0, 1, num_rv_samples); repmat(ref, 1, num_rv_samples)];
end

var_data_cv = zeros(length(num_estimator_samples), 1);
var_data_acv = zeros(length(num_estimator_samples), 1);
for i=1:length(num_estimator_samples)
  est_samples = num_estimator_samples(i);
  fprintf("\n*** Est Samples: %d ***\n", est_samples);
  temp_data_cv = zeros(est_samples, 1);
  temp_data_acv = zeros(est_samples, 1);
  for j=1:est_samples
    temp_data_cv(j) = cv.est(x0_rv_ext(:, :, j), n_cv, Uopt_h_num, false);
    temp_data_acv(j) = acv.est(x0_rv_ext(:, :, j), n_acv, m_acv, Uopt_h_num, false, 'stat', 'stat');
  end
  var_data_cv(i) = var(temp_data_cv);
  var_data_acv(i) = var(temp_data_acv);
end

%% plot
var_cv_an = cv.variance(n_cv, Uopt_h_num);
var_acv_an = acv.variance(n_acv, m_acv, Uopt_h_num);
se_cv = (var_data_cv - var_cv_an).^2;
se_acv = (var_data_acv - var_acv_an).^2;
figure;
hold on;
semilogy(num_estimator_samples, se_cv, 'g', 'LineWidth', 2, 'DisplayName', 'CV Var squared error');
semilogy(num_estimator_samples, se_acv, 'r', 'LineWidth', 2, 'DisplayName', 'ACV Var squared error');
xlabel("Estimator Samples");
ylabel("Variance");
title("Variance difference compared to analytical ($\alpha=$stat)", "Interpreter", "latex");
legend show;
grid on;

%% G Convergence and variance with various optimizers and sample sizes
num_rv_samples = [10 100 500];
num_rv_samples_actual = zeros(length(num_rv_samples), 4); % 4 bc we have [MC, CV, n_ACV ,m_ACV]
num_estimator_samples = 50;

u0_num = repelem(data.lqrsol{2}.Uopt, 10, 1); % warm start
max_iters = 30;
tol = 1e-12;

lf_hf_cost_ratio = 0.2; % TODO get this experimentally for each num_samples
acv_exp_sample_factor = 1; % factor to increase samples for ACV Expectation, m

mc = Mc(x0_ext_mean, x0_ext_cov, data.lqrsol{1});
cv = Cv(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2});
acv = Acv(x0_ext_mean, x0_ext_cov, data.lqrsol{1}, data.lqrsol{2});

data.hf_cost = zeros(max_iters, num_estimator_samples, length(num_rv_samples));
data.hf_u = zeros(length(u0_num), max_iters, num_estimator_samples, length(num_rv_samples));
data.cv_cost = zeros(max_iters, num_estimator_samples, length(num_rv_samples));
data.cv_u = zeros(length(u0_num), max_iters, num_estimator_samples, length(num_rv_samples));
data.acv_cost = zeros(max_iters, num_estimator_samples, length(num_rv_samples));
data.acv_u = zeros(length(u0_num), max_iters, num_estimator_samples, length(num_rv_samples));

%%
for num_samples=num_rv_samples
  fprintf("\n*** RV Samples: %d ***\n", num_samples);
  wait_bar = waitbar(0, "Running estimators with " + num_samples + " samples");
  n_mc = num_samples;
  n_cv = int32(n_mc / (1 + lf_hf_cost_ratio));
  n_acv = round(num_samples / (1 + lf_hf_cost_ratio + acv_exp_sample_factor*lf_hf_cost_ratio));
  m_acv = round(acv_exp_sample_factor * n_acv);
  num_rv_samples_actual(num_rv_samples == num_samples, :) = [n_mc, n_cv, n_acv, m_acv];
  
  n_total = round(max(n_cv+n_acv+m_acv, num_samples)); % account for weird splits like 0.999
  
  for i=1:num_estimator_samples
    waitbar(i/num_estimator_samples, wait_bar);
    % get the RV samples
    x0_rv = mvnrnd(x0_mean, x0_cov, n_total)';
    x0_rv_ext = [x0_rv; repmat(u0, 1, n_total); repmat(ref, 1, n_total)];
    
    % MC with HF
    [costs, Us, ~] = mc.opt(u0_num, max_iters, tol, x0_rv_ext, n_mc);
    data.hf_cost(:, i, num_rv_samples == num_samples) = costs;
    data.hf_u(:, :, i, num_rv_samples == num_samples) = Us;
    
    % CV (true = lf at max corr)
    [costs, Us, ~] = cv.opt(u0_num, max_iters, tol, x0_rv_ext, n_cv, true);
    data.cv_cost(:, i, num_rv_samples == num_samples) = costs;
    data.cv_u(:, :, i, num_rv_samples == num_samples) = Us;
    
    % ACV (true = lf at max corr)
    [costs, Us, ~] = acv.opt(u0_num, max_iters, tol, x0_rv_ext, n_acv, m_acv, true, '-1', 'stat');
    data.acv_cost(:, i, num_rv_samples == num_samples) = costs;
    data.acv_u(:, :, i, num_rv_samples == num_samples) = Us;
  end
  close(wait_bar);
end


%% plot percentiles
percentiles = [25 75];
markers = ["-", "--"];
data.hf_prct = prctile(data.hf_cost, percentiles, 2);
data.cv_prct = prctile(data.cv_cost, percentiles, 2);
data.acv_prct = prctile(data.acv_cost, percentiles, 2);
fig = figure;
fig.Position(3:4) = [1500 500];
sgtitle("Cost Percentiles vs iterations");
for i=1:length(num_rv_samples)
  subplot(1, length(num_rv_samples), i);
  hold on;
  for j=1:length(percentiles)
    semilogy(5:max_iters, data.hf_prct(5:end, j, i), 'r', 'LineWidth', 1, 'LineStyle', markers(j), 'DisplayName', "MC HF "+percentiles(j)+"%");
    semilogy(5:max_iters, data.cv_prct(5:end, j, i), 'g', 'LineWidth', 1, 'LineStyle', markers(j), 'DisplayName', "CV "+percentiles(j)+"%");
    semilogy(5:max_iters, data.acv_prct(5:end, j, i), 'b', 'LineWidth', 1, 'LineStyle', markers(j), 'DisplayName', "ACV "+percentiles(j)+"%");
  end
  
  title("Samples: MC="+num_rv_samples_actual(i, 1)+", CV n="+num_rv_samples_actual(i, 2)+", ACV n="+num_rv_samples_actual(i, 3)+", m="+num_rv_samples_actual(i, 4));
  xlabel("Iteration");
  ylabel("Cost");
  legend show;
  % ylim([1e2 2e2]);
  grid on;
end

%% plot variance and box plots
% 2nd dimension is estimator samples, squeeze to get rid of extra dimension
data.hf_var_st = squeeze(var(data.hf_cost, 0, 2));
data.cv_var_st = squeeze(var(data.cv_cost, 0, 2));
data.acv_var_st = squeeze(var(data.acv_cost, 0, 2));
fig = figure;
fig.Position(3:4) = [1500 500];
sgtitle("Variance of cost vs iterations");
for i=1:length(num_rv_samples)
  subplot(1, length(num_rv_samples), i);
  hold on;
  
  % yyaxis left;
  semilogy(1:max_iters, data.hf_var_st(:, i), 'r', 'LineWidth', 1, 'DisplayName', 'MC HF');
  semilogy(1:max_iters, data.cv_var_st(:, i), 'g', 'LineWidth', 1, 'DisplayName', 'CV');
  semilogy(1:max_iters, data.acv_var_st(:, i), 'b', 'LineWidth', 1, 'DisplayName', 'ACV');
  ylabel("Variance");
  
  % yyaxis right;
  % boxplot(data.hf_cost(5:max_iters, :, i)', 'Positions', 5:max_iters, 'Colors', 'r', 'Widths', 0.5, 'Whisker', 1);
  % boxplot(data.cv_fix_cost(5:max_iters, :, i)', 'Positions', 5:max_iters, 'Colors', 'g', 'Widths', 1, 'Whisker', 1);
  % boxplot(data.acv_fix_cost(5:max_iters, :, i)', 'Positions', 5:max_iters, 'Colors', 'b', 'Widths', 0.5, 'Whisker', 1);
  % ylabel("Cost/Objective");
  
  ax = gca;
  ax.YAxis(1).Scale ="log";
  % ax.YAxis(2).Scale ="log";
  % ax.XAxis.TickValues = 1:max_iters;
  title("Samples: MC="+num_rv_samples_actual(i, 1)+", CV n="+num_rv_samples_actual(i, 2)+", ACV n="+num_rv_samples_actual(i, 3)+", m="+num_rv_samples_actual(i, 4));
  xlabel("Iteration");
  legend show;
  % ylim([5e0 2e5]);
  grid on;
end

%% plot convergence in solution
for i=1:length(num_rv_samples)
  plotConvergence(data.lqrsol{1}.Uopt, squeeze(data.hf_u(:, :, :, i)), squeeze(data.cv_u(:, :, :, i)), squeeze(data.acv_u(:, :, :, i)), ...
    ["MC HF", "CV", "ACV"], "Convergence distance comparison with " + num_rv_samples(i) + " HF samples");
end

%% save all data
save("artefacts/data.mat");

function xdot = msd(~, x, u, A, B)
xdot = A * x + B * u;
end

function x_history = fminuncWHistory(fun, x0)
% this is a wrapper around fminunc that returns the history of the optimization
x_history = [];
options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'OutputFcn', @outfun);
fminunc(fun, x0, options);

  function stop = outfun(x, ~, state)
    % this is for capturing intermediate values of the optimization
    stop = false;
    if isequal(state, 'iter')
      x_history = [x_history x];
    end
  end
end

function x = simCL(Tsim, Tfid, Ts, x0, u0, U, A, B, dyn)
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

function perturbation = perturbMaxGrad(lqrsol, perturb_dir_samples, perturb_dir_mag, perturb_range)
rand_grads = zeros(length(lqrsol.Uopt), perturb_dir_samples);
% generate uniform random samples in range [-1, 1]
perturb_dirs = - 1 + 2.*rand(length(lqrsol.Uopt), perturb_dir_samples);
% normalize dirs so they are unit vectors
perturb_dirs = perturb_dirs ./ sqrt(sum(perturb_dirs.^2, 1));
% scale the vector to have a magnitude of perturb_dir_mag
perturb_dirs = perturb_dirs * perturb_dir_mag;
% check the gradient in random directions
for i = 1:perturb_dir_samples
  rand_grads(:, i) = St.LQRGrad(lqrsol.x0_ext, lqrsol, lqrsol.Uopt + perturb_dirs(:, i));
end
% pick max gradient direction
[~, max_grad_idx] = max(sum(rand_grads.^2, 1));
perturb_dir_max = perturb_dirs(:, max_grad_idx)./ perturb_dir_mag; % length 1 so can scale later

% perturb the solution in the max gradient direction
perturbation = perturb_dir_max * perturb_range;
end

function plotMultifidCost(perturb_range, cost_hf, cost_lf, corr_st, corr_an, title_str, costs_str, xaxis_str)
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

function plotCorr2d(corr, title_str, sols_str, save_str)
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

function visSols(uopt_hf, uopt_lf, times, title_str, zaxis, zaxis_str)
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
sld.ValueChangedFcn = @(sld, event) updVisSol(sld, plot_hf, plot_lf, plot_hf_raw, plot_lf_raw, uopt_hf, uopt_lf);
end

function updVisSol(sld, plot_hf, plot_lf, plot_hf_raw, plot_lf_raw, uopt_hf, uopt_lf)
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

function plotConvergence(an_sol, U_num_hf, U_num_cv, U_num_cv_fix, labels, title_str)
dist_hf = squeeze(vecnorm(an_sol - U_num_hf, 2, 1));
dist_cv = squeeze(vecnorm(an_sol - U_num_cv, 2, 1));
dist_cv_fix = squeeze(vecnorm(an_sol - U_num_cv_fix, 2, 1));
if size(dist_hf, 1) == 1
  dist_hf_mean = dist_hf;
  dist_cv_mean = dist_cv;
  dist_cv_fix_mean = dist_cv_fix;
else
  dist_hf_mean = mean(dist_hf, 2);
  dist_cv_mean = mean(dist_cv, 2);
  dist_cv_fix_mean = mean(dist_cv_fix, 2);
end
figure;
% handle alpha to be min 0.1 and max 1 based on how many lines there are
semilogy(dist_hf, 'LineWidth', 1, 'Color', [1 0 0, max(0.1, 1/size(dist_hf, 1))]);
hold on;
% plot means with thicker lines
semilogy(dist_cv, 'LineWidth', 1, 'Color', [0 1 0, max(0.1, 1/size(dist_cv, 1))]);
semilogy(dist_cv_fix, 'LineWidth', 1, 'Color', [0 0 1, max(0.1, 1/size(dist_cv_fix, 1))]);
l{1} = semilogy(dist_hf_mean, 'r', "Marker", "o", 'LineWidth', 1.5);
l{2} = semilogy(dist_cv_mean, 'g', "Marker", "o", 'LineWidth', 1.5);
l{3} = semilogy(dist_cv_fix_mean, 'b', "Marker", "o", 'LineWidth', 1.5);
title(title_str);
xlabel("Iteration");
ylabel("Distance to analytical solution");
legend([l{:}], labels);
grid on;
end
