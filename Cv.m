classdef Cv < handle
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
		function obj = Cv(x0_mean, x0_cov, lqrsol_hf, lqrsol_lf, l_h_cost_ratio)
			obj.x0_mean = x0_mean;
			obj.x0_cov = x0_cov;
			obj.lqrsol_hf = lqrsol_hf;
			obj.lqrsol_lf = lqrsol_lf;
			obj.l_h_cost_ratio = l_h_cost_ratio;
		end
		
		function cost = est(obj, x0_rv_ext, n, u, use_best_U_lf)
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
			var_l = St.LQRVar(obj.x0_mean, obj.x0_cov, obj.lqrsol_lf, u_hla); % analytical
			cov_hl = St.LQRCov(obj.x0_mean, obj.x0_cov, obj.lqrsol_hf, obj.lqrsol_lf, u, u_hla); % analytical
			exp_l = St.LQRExp(obj.x0_mean, obj.x0_cov, obj.lqrsol_lf, u_hla);
			alpha = -cov_hl / var_l;
			cost = cost_hf + alpha * (cost_lf - exp_l);
		end
		
		function [costs, Us, U_hlas] = opt(obj, u0, max_iters, tol, x0_rv_ext, n, use_best_U_lf, use_sgd)
			costs = zeros(max_iters, 1);
			obj.Us = zeros(size(u0, 1), max_iters);
			obj.U_hlas = zeros(size(u0, 1)/10, max_iters);
			obj.costs_lf = zeros(n, max_iters); % we use only n samples out of all in x0_rv
			obj.idx = 1;
			
			if max_iters <  0 || tol < 0
				options = optimoptions('fminunc', 'SpecifyObjectiveGradient', false, 'OutputFcn', @OutFn);
				% options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'OutputFcn', @OutFn);
			else
				options = optimoptions('fminunc', 'SpecifyObjectiveGradient', false, 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', tol, 'StepTolerance', tol);
				% options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', tol, 'StepTolerance', tol, 'Display', 'iter-detailed');
			end
			f = @(u) obj.est(x0_rv_ext, n, u, use_best_U_lf);
			
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
				% function [f, g] = fwGrad(x0_rv_ext, n, u, use_best_U_lf)
				% 	f = obj.est(x0_rv_ext, n, u, use_best_U_lf);
				% 	g = mean(St.LQRGrad(x0_rv_ext(:, 1:n), obj.lqrsol_hf, u), 2);
				% 	% fprintf('grad: ');
				% 	% disp(g);
				% end
				% f = @(u) fwGrad(x0_rv_ext, n, u, use_best_U_lf);
				% grad_opts = optimoptions("fminunc", FiniteDifferenceType="central");
				% checkGradients(f, u0, grad_opts, Display="on");
				fminunc(f, u0, options);
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
		
		function var = variance(obj, n, u)
			var_h = St.LQRVar(obj.x0_mean, obj.x0_cov, obj.lqrsol_hf, u);
			u_lf = St.DownsampleAvg(u, 10);
			corr_hl = St.LQRCorr(obj.x0_mean, obj.x0_cov, obj.lqrsol_hf, obj.lqrsol_lf, u, u_lf);
			var = var_h/n * (1 - corr_hl^2);
		end
		
		function n_cv = getEqCostSamples(obj, n_mc)
			n_cv = round(n_mc / (1 + obj.l_h_cost_ratio));
		end
		
	end
	
end