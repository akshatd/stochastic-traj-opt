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
			cost = mean(St.LQRCost(x0_rv_ext(:, 1:n), obj.lqrsol, u));
		end
		
		function [costs, Us, U_hlas] = opt(obj, u0, max_iters, tol, x0_rv_ext, n)
			costs = zeros(max_iters, 1);
			obj.Us = zeros(size(u0, 1), max_iters);
			obj.U_hlas = zeros(size(u0, 1)/10, max_iters);
			obj.idx = 1;
			
			if max_iters <  0 || tol < 0
				options = optimoptions('fminunc', 'OutputFcn', @OutFn);
			else
				options = optimoptions('fminunc', 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', tol);
			end
			fminunc(@(u) obj.est(x0_rv_ext, n, u), u0, options);
			
			% trim to match iters
			obj.idx = obj.idx-1; % remove last iteration that stopped it
			if max_iters < 0
				iters = obj.idx;
			else
				iters = min(obj.idx, max_iters);
			end
			costs = costs(1:iters);
			obj.Us = obj.Us(:, 1:iters);
			obj.U_hlas = Est.DownsampleAvg(obj.Us, 10);
			
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