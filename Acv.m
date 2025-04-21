classdef Acv < handle
	properties
		Us = [];
		U_hlas = [];
		costs_lf = [];
		idx = 1;
		x0_mean = [];
		x0_cov = [];
		lqrsol_hf = {};
		lqrsol_lf = {};
		l_h_cost_ratio = 0.045;
	end
	
	methods
		function obj = Acv(x0_mean, x0_cov, lqrsol_hf, lqrsol_lf, l_h_cost_ratio)
			obj.x0_mean = x0_mean;
			obj.x0_cov = x0_cov;
			obj.lqrsol_hf = lqrsol_hf;
			obj.lqrsol_lf = lqrsol_lf;
			obj.l_h_cost_ratio = l_h_cost_ratio;
		end
		
		function cost = est(obj, x0_rv_ext, n, m, u, use_best_U_lf, a_type, e_type)
			% ONLY set use_best_U_lf if insider optimizer and idx is set outside
			cost_hf_all = St.LQRObj(x0_rv_ext(:, 1:n), obj.lqrsol_hf, u);
			u_hla = St.DownsampleAvg(u, 10);
			cost_lf_all = St.LQRObj(x0_rv_ext(:, 1:n), obj.lqrsol_lf, u_hla);
			if use_best_U_lf
				obj.Us(:, obj.idx) = u;
				obj.costs_lf(:, obj.idx) = cost_lf_all;
				% TODO: Us should be in rows to prevent transpose
				corrs = St.CorrMulti2D(cost_hf_all', obj.costs_lf(:, 1:obj.idx)');
				[~, best_idx] = max(corrs);
				u_hla = St.DownsampleAvg(obj.Us(:, best_idx), 10);
				cost_lf_all = obj.costs_lf(:, best_idx);
			end
			obj.U_hlas(:, obj.idx) = u_hla; % save bc used for plotting
			cost_hf = mean(cost_hf_all);
			cost_lf = mean(cost_lf_all);
			if strcmp(a_type, 'anly')
				var_l = St.LQRVar(obj.x0_mean, obj.x0_cov, obj.lqrsol_lf, u_hla); % analytical
				cov_hl = St.LQRCov(obj.x0_mean, obj.x0_cov, obj.lqrsol_hf, obj.lqrsol_lf, u, u_hla); % analytical
				alpha = -m/(m+n) * cov_hl / var_l;
			elseif strcmp(a_type, 'stat')
				var_l = var(cost_lf_all);
				cov_hl = cov(cost_hf_all, cost_lf_all);
				cov_hl = cov_hl(1, 2); % only off diagonal element
				% alpha = -m/(m+n) * cov_hl / var_l; % based off notes
				alpha = -cov_hl / var_l; % based off paper
			elseif strcmp(a_type, '-1')
				alpha = -1;
				% alpha = -5.8;
				% alpha = -5.2;
			end
			if strcmp(e_type, 'anly')
				exp_l = St.LQRExp(obj.x0_mean, obj.x0_cov, obj.lqrsol_lf, u_hla);
			elseif strcmp(e_type, 'stat')
				exp_l = mean(St.LQRObj(x0_rv_ext(:, n+1:n+m), obj.lqrsol_lf, u_hla));
			elseif strcmp(e_type, 'share')
				exp_l = St.LQRObj(x0_rv_ext(:, n+1:n+m), obj.lqrsol_lf, u_hla);
				exp_l = mean([exp_l; cost_lf_all]);
			end
			cost = cost_hf + alpha * (cost_lf - exp_l);
		end
		
		function cost = estPrecalc(obj, x0_rv_ext, n, m, u, use_best_U_lf, a_type, e_type, x0_term_hf, x0_term_lf, x0_term_exp)
			% ONLY set use_best_U_lf if insider optimizer and idx is set outside
			% cost_hf_all = St.LQRObj(x0_rv_ext(:, 1:n), obj.lqrsol_hf, u);
			cost_hf_all = St.LQRObj_precalc(x0_rv_ext(:, 1:n), obj.lqrsol_hf, u, x0_term_hf);
			u_hla = St.DownsampleAvg(u, 10);
			% cost_lf_all = St.LQRObj(x0_rv_ext(:, 1:n), obj.lqrsol_lf, u_hla);
			cost_lf_all = St.LQRObj_precalc(x0_rv_ext(:, 1:n), obj.lqrsol_lf, u_hla, x0_term_lf);
			if use_best_U_lf
				obj.Us(:, obj.idx) = u;
				obj.costs_lf(:, obj.idx) = cost_lf_all;
				% TODO: Us should be in rows to prevent transpose
				corrs = St.CorrMulti2D(cost_hf_all', obj.costs_lf(:, 1:obj.idx)');
				[~, best_idx] = max(corrs);
				u_hla = St.DownsampleAvg(obj.Us(:, best_idx), 10);
				cost_lf_all = obj.costs_lf(:, best_idx);
			end
			obj.U_hlas(:, obj.idx) = u_hla; % save bc used for plotting
			cost_hf = mean(cost_hf_all);
			cost_lf = mean(cost_lf_all);
			if strcmp(a_type, 'anly')
				var_l = St.LQRVar(obj.x0_mean, obj.x0_cov, obj.lqrsol_lf, u_hla); % analytical
				cov_hl = St.LQRCov(obj.x0_mean, obj.x0_cov, obj.lqrsol_hf, obj.lqrsol_lf, u, u_hla); % analytical
				alpha = -m/(m+n) * cov_hl / var_l;
			elseif strcmp(a_type, 'stat')
				var_l = var(cost_lf_all);
				cov_hl = cov(cost_hf_all, cost_lf_all);
				cov_hl = cov_hl(1, 2); % only off diagonal element
				% alpha = -m/(m+n) * cov_hl / var_l; % based off notes
				alpha = -cov_hl / var_l; % based off paper
			elseif strcmp(a_type, '-1')
				alpha = -1;
				% alpha = -5.8;
				% alpha = -5.2;
			end
			if strcmp(e_type, 'anly')
				exp_l = St.LQRExp(obj.x0_mean, obj.x0_cov, obj.lqrsol_lf, u_hla);
			elseif strcmp(e_type, 'stat')
				% exp_l = mean(St.LQRObj(x0_rv_ext(:, n+1:n+m), obj.lqrsol_lf, u_hla));
				exp_l = mean(St.LQRObj_precalc(x0_rv_ext(:, n+1:n+m), obj.lqrsol_lf, u_hla, x0_term_exp));
			elseif strcmp(e_type, 'share')
				% exp_l = St.LQRObj(x0_rv_ext(:, n+1:n+m), obj.lqrsol_lf, u_hla);
				exp_l = St.LQRObj_precalc(x0_rv_ext(:, n+1:n+m), obj.lqrsol_lf, u_hla, x0_term_exp);
				exp_l = mean([exp_l; cost_lf_all]);
			end
			cost = cost_hf + alpha * (cost_lf - exp_l);
		end
		
		
		function [costs, Us, U_hlas] = opt(obj, u0, max_iters, tol, x0_rv_ext, n, m, use_best_U_lf, a_type, e_type, use_sgd, U_bounds_A, U_bounds_b)
			costs = zeros(max_iters, 1);
			obj.Us = zeros(size(u0, 1), max_iters);
			obj.U_hlas = zeros(size(u0, 1)/10, max_iters);
			obj.costs_lf = zeros(n, max_iters); % we use only n samples out of all in x0_rv
			obj.idx = 1;
			
			if max_iters <  0 || tol < 0
				options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'OutputFcn', @OutFn);
				% options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'OutputFcn', @OutFn);
			else
				options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', tol, 'StepTolerance', tol);
				% options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', tol, 'StepTolerance', tol, 'Display', 'iter-detailed');
			end
			% normal f
			% f = @(u) obj.est(x0_rv_ext, n, m, u, use_best_U_lf, a_type, e_type);
			
			% f with precalculated x0 term
			x0_term_hf = St.LQRObj_x0term(x0_rv_ext(:, 1:n), obj.lqrsol_hf);
			x0_term_lf = St.LQRObj_x0term(x0_rv_ext(:, 1:n), obj.lqrsol_lf);
			x0_term_exp = St.LQRObj_x0term(x0_rv_ext(:, n+1:n+m), obj.lqrsol_lf);
			f = @(u) obj.estPrecalc(x0_rv_ext, n, m, u, use_best_U_lf, a_type, e_type, x0_term_hf, x0_term_lf, x0_term_exp);
			% SGD algorithm
			if use_sgd
				% alpha_start = 7e-9; % learning rate
				alpha_start = 2e-4; % normie learning rate
				u = u0; % initial control input
				for i=1:max_iters
					% Store the cost
					costs(i) = f(u);
					obj.Us(:, obj.idx) = u;
					obj.U_hlas(:, obj.idx) = St.DownsampleAvg(u, 10);
					% Compute the gradient
					x0_idxs = randperm(n);
					x0_samples = x0_rv_ext(:, x0_idxs(1:round(n/2)));
					grad = mean(St.LQRGrad(x0_samples, obj.lqrsol_hf, u), 2);
					% disp(norm(grad));
					% Update the control input
					% alpha = alpha_start*norm(grad);
					alpha = alpha_start;
					u = u - alpha * grad;
					obj.idx = obj.idx + 1;
					% Check for convergence
					if i > 1 && abs(costs(i) - costs(i-1)) < tol
						break;
					end
				end
			else
				% function [f, g] = fwGrad(x0_rv_ext, n, m, u, use_best_U_lf, a_type, e_type)
				% 	f = obj.est(x0_rv_ext, n, m, u, use_best_U_lf, a_type, e_type);
				% 	g = mean(St.LQRGrad(x0_rv_ext(:, 1:n), obj.lqrsol_hf, u), 2);
				% 	% fprintf('grad: ');
				% 	% disp(g);
				% end
				% f = @(u) fwGrad(x0_rv_ext, n, m, u, use_best_U_lf, a_type, e_type);
				% grad_opts = optimoptions("fminunc", FiniteDifferenceType="central");
				% checkGradients(f, u0, grad_opts, Display="on");
				% fminunc(f, u0, options);
				fmincon(f, u0, U_bounds_A, U_bounds_b, [], [], [], [], [], options);
			end
			
			% trim to match iters
			obj.idx = obj.idx-1; % remove last iteration that stopped it
			if max_iters < 0
				iters = obj.idx;
			else
				iters = min(obj.idx, max_iters);
			end
			costs = costs(1:iters);
			obj.Us = obj.Us(:, 1:iters);
			obj.U_hlas = obj.U_hlas(:, 1:iters);
			
			function stop = OutFn(x, optimValues, state)
				stop = false;
				if isequal(state, 'iter')
					% obj.idx starts at 1, so increment after assignment
					obj.Us(:, obj.idx) = x;
					costs(obj.idx) = optimValues.fval;
					obj.idx = obj.idx+1;
				end
			end
			Us = obj.Us;
			U_hlas = obj.U_hlas;
		end
		
		function var = variance(obj, n, m, u)
			var_h = St.LQRVar(obj.x0_mean, obj.x0_cov, obj.lqrsol_hf, u);
			u_lf = St.DownsampleAvg(u, 10);
			corr_hl = St.LQRCorr(obj.x0_mean, obj.x0_cov, obj.lqrsol_hf, obj.lqrsol_lf, u, u_lf);
			% var = var_h/n * (1 - m/(m+n) * corr_hl^2); % for MLMC ACV, notes
			r1 = m/n;
			var = var_h/n * (1 - (r1-1)/r1 * corr_hl^2); % for MFMC ACV, paper
		end
		
		function [n_acv, m_acv] = getEqCostSamples(obj, n_mc, mn_ratio)
			n_acv = round(n_mc / (1 + obj.l_h_cost_ratio + mn_ratio*obj.l_h_cost_ratio));
			m_acv = round(mn_ratio * n_acv);
		end
		
		% following functions are only for optimization
		function var = varianceEqCost(obj, n_cv, mn_ratio, u)
			[n, m] = obj.getEqCostSamplesRaw(n_cv, mn_ratio);
			var = obj.variance(n, m, u);
		end
		
		function [n_acv, m_acv] = getEqCostSamplesRaw(obj, n_mc, mn_ratio)
			n_acv = n_mc / (1 + obj.l_h_cost_ratio + mn_ratio*obj.l_h_cost_ratio);
			m_acv = mn_ratio * n_acv;
		end
		
	end
end