classdef Mc < handle
	properties
		Us = [];
		U_hlas = [];
		idx = 1;
		x0_mean = [];
		x0_cov = [];
		lqrsol = {};
	end
	
	methods
		function obj = Mc(x0_mean, x0_cov, lqrsol)
			obj.x0_mean = x0_mean;
			obj.x0_cov = x0_cov;
			obj.lqrsol = lqrsol;
		end
		
		function cost = est(obj, x0_rv_ext, n, u)
			cost = mean(St.LQRObj(x0_rv_ext(:, 1:n), obj.lqrsol, u));
		end
		
		function cost = estPrecalc(obj, x0_rv_ext, n, u, x0_term)
			cost = mean(St.LQRObj_precalc(x0_rv_ext(:, 1:n), obj.lqrsol, u, x0_term));
		end
		
		function [costs, Us, U_hlas] = opt(obj, u0, max_iters, tol, x0_rv_ext, n, use_sgd, U_bounds_A, U_bounds_b)
			costs = zeros(max_iters, 1);
			obj.Us = zeros(size(u0, 1), max_iters);
			obj.U_hlas = zeros(size(u0, 1)/10, max_iters);
			obj.idx = 1;
			
			if max_iters <  0 || tol < 0
				options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'OutputFcn', @OutFn);
				% options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'OutputFcn', @OutFn);
			else
				options = optimoptions('fmincon', 'SpecifyObjectiveGradient', false, 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', tol, 'StepTolerance', tol);
				% options = optimoptions('fminunc', 'SpecifyObjectiveGradient', true, 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', tol, 'StepTolerance', tol, 'Display', 'iter-detailed');
			end
			
			% normal f
			% f = @(u) obj.est(x0_rv_ext, n, u);
			
			% f with precalculated x0 term
			x0_term = St.LQRObj_x0term(x0_rv_ext(:, 1:n), obj.lqrsol);
			f = @(u) obj.estPrecalc(x0_rv_ext, n, u, x0_term);
			
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
					grad = mean(St.LQRGrad(x0_samples, obj.lqrsol, u), 2);
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
				% function [f, g] = fwGrad(x0_rv_ext, n, u)
				% 	f = obj.est(x0_rv_ext, n, u);
				% 	g = mean(St.LQRGrad(x0_rv_ext(:, 1:n), obj.lqrsol, u), 2);
				% 	% fprintf('grad: ');
				% 	% disp(g);
				% end
				% f = @(u) fwGrad(x0_rv_ext, n, u);
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
			obj.U_hlas = St.DownsampleAvg(obj.Us, 10);
			
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
		
		function var = var(n, u)
			var_h = St.LQRVar(obj.x0_mean, obj.x0_cov, obj.lqrsol, u);
			var = var_h/n;
		end
		
		
	end
	
end