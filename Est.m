% class for all the estimators Control variate, approximate control variate, monte carlo etc

classdef Est
	methods(Static)
		function [costs, Us] = CvOpt(u0, max_iters, x0_rv_ext, x0_mean, x0_cov, lqrsol_hf, lqrsol_lf, n, use_best_U_lf)
			costs = zeros(max_iters, 1);
			Us = zeros(size(u0, 1), max_iters);
			costs_lf = zeros(n, max_iters); % we use only n samples out of all in x0_rv
			curr_iter = 1;
			
			options = optimoptions('fminunc', 'OutputFcn', @OutFn, 'MaxIter', max_iters, 'OptimalityTolerance', 1e-12);
			fminunc(@(u) CvEst(x0_rv_ext, x0_mean, x0_cov, lqrsol_hf, lqrsol_lf, n, u, use_best_U_lf), u0, options);
			
			% trim bc it somehow runs for more iterations
			costs = costs(1:max_iters);
			Us = Us(:, 1:max_iters);
			function cost = CvEst(x0_rv_ext, x0_mean, x0_cov, lqrsol_hf, lqrsol_lf, n, u, use_best_U_lf)
				% MC estimator for hf
				cost_hf_all = St.LQRCost(x0_rv_ext(:, 1:n), lqrsol_hf, u);
				cost_hf = mean(cost_hf_all);
				% reuse same samples for LF
				u_hla = Est.DownsampleAvg(u, 10);
				cost_lf_all = St.LQRCost(x0_rv_ext(:, 1:n), lqrsol_lf, u_hla);
				if use_best_U_lf
                    Us(:, curr_iter) = u; % temporary init so it can be used here
                    costs_lf(:, curr_iter) = cost_lf_all;

					% TODO: Us should be in rows to prevent transpose
					corrs = St.CorrMulti2D(cost_hf_all', costs_lf(:, 1:curr_iter)');
					[~, best_idx] = max(corrs);
					% fprintf('It: %d, Best idx: %d, Corr: %f\n', curr_iter, best_idx, corrs(best_idx));
					u_hla = Est.DownsampleAvg(Us(:, best_idx), 10);
					cost_lf_all = costs_lf(:, best_idx);
				end
				cost_lf = mean(cost_lf_all);
				% use actual mean, cov for expectation
				exp_l = St.LQRExp(x0_mean, x0_cov, lqrsol_lf, u_hla);
				var_l = St.LQRVar(x0_mean, x0_cov, lqrsol_lf, u_hla);
				% calculate optimal alpha
				cov_hl = St.LQRCov(x0_mean, x0_cov, lqrsol_hf, lqrsol_lf, u, u_hla);
				alpha = -cov_hl / var_l;
				cost = cost_hf + alpha * (cost_lf - exp_l);
			end
			
			function stop = OutFn(x, optimValues, state)
				stop = false;
				if isequal(state, 'iter')
                    % curr_iter starts at 1, so increment after assignment
                    Us(:, curr_iter) = x;
					costs(curr_iter) = optimValues.fval;
					curr_iter = curr_iter+1;
				end
			end
		end
		
		function U_lf = DownsampleAvg(U_hf, avg_window, repeat)
			arguments
				U_hf (:, :) double
				avg_window int32
				repeat logical = false
			end
			U_hf_len = size(U_hf, 1);
			U_lf = zeros(U_hf_len/avg_window, size(U_hf, 2));
			skips = U_hf_len/avg_window;
			for i = 1:skips
				row_start = (i-1) * avg_window + 1;
				row_end = i*avg_window;
				U_lf(i, :) = mean(U_hf(row_start:row_end, :), 1);
			end
			if repeat
				U_lf = repelem(U_lf, avg_window, 1);
			end
		end
		
	end
end