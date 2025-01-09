% class that contains functions relevant to calculating statistics

classdef St
	methods(Static) % so that we don't need to create an object to use these functions
		function cost = LQRCost(x0, U, lqrsol)
			Q = lqrsol.Q; S = lqrsol.S; M = lqrsol.M; Qbar = lqrsol.Qbar; Rbar = lqrsol.Rbar;
			% handles x0 having multiple samples with diag
			cost = U'*(S'*Qbar*S + Rbar)*U + 2*x0'*M'*Qbar*S*U + diag(x0'*(M'*Qbar*M + Q)*x0);
		end
		
		function grad = LQRGrad(x0, U, lqrsol)
			S = lqrsol.S; M = lqrsol.M; Qbar = lqrsol.Qbar; Rbar = lqrsol.Rbar;
			% 2HU + 2q
			H = S'*Qbar*S + Rbar;
			q = x0'*M'*Qbar*S;
			grad = 2 * H * U + 2 * q';
		end
		
		function [f, g] = LQRCostwGrad(x0, U, lqrsol)
			f = St.LQRCost(x0, U, lqrsol);
			g = St.LQRGrad(x0, U, lqrsol);
		end
		
		function exp = LQRExp(x0_mean, x0_cov, U, Q, S, M, Qbar, Rbar)
			K = U'*(S'*Qbar*S + Rbar)*U;
			L = 2 * M' * Qbar * S * U;
			N = M' * Qbar * M + Q;
			exp = K + x0_mean' * L + x0_mean' * N * x0_mean + trace(N * x0_cov);
		end
		
		function var = LQRVar(x0_mean, x0_cov, U, Q, S, M, Qbar)
			L = 2 * M' * Qbar * S * U;
			N = M' * Qbar * M + Q;
			var = L' * x0_cov * L + 2 * trace(N * x0_cov * N * x0_cov) + 4 * (x0_mean' * N + L') * x0_cov * N * x0_mean;
		end
		
		function lqr_cov = LQRCov(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, Uopt_hf, Uopt_lf)
			% K1 = Uopt_hf' * (lqrsol_hf.S' * lqrsol_hf.Qbar * lqrsol_hf.S + lqrsol_hf.Rbar) * Uopt_hf;
			L1 = 2 * lqrsol_hf.M' * lqrsol_hf.Qbar * lqrsol_hf.S * Uopt_hf;
			N1 = lqrsol_hf.M' * lqrsol_hf.Qbar * lqrsol_hf.M + lqrsol_hf.Q;
			% K2 = Uopt_lf' * (lqrsol_lf.S' * lqrsol_lf.Qbar * lqrsol_lf.S + lqrsol_lf.Rbar) * Uopt_lf;
			L2 = 2 * lqrsol_lf.M' * lqrsol_lf.Qbar * lqrsol_lf.S * Uopt_lf;
			N2 = lqrsol_lf.M' * lqrsol_lf.Qbar * lqrsol_lf.M + lqrsol_lf.Q;
			x0_rv_ext = [x0_rv; repmat(u0, 1, size(x0_rv, 2)); repmat(ref, 1, size(x0_rv, 2))];
			x_mean = mean(x0_rv_ext, 2);
			x_cov = cov(x0_rv_ext');
			lqr_cov = L1'*x_cov*L2 + 2*x_mean'*N2*x_cov*L1 + 2*x_mean'*N1*x_cov*L2 + 2*trace(N1*x_cov*N2*x_cov) + 4*x_mean'*N1*x_cov*N2*x_mean;
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
				cov_an = St.LQRCov(x0_rv, u0, ref, lqrsol_hf, lqrsol_lf, U_hf(:, i), U_lf(:, i)); % TODO make indices consistent
				var_J1 = St.LQRVar(x_mean, x_cov, U_hf(:, i), lqrsol_hf.Q, lqrsol_hf.S, lqrsol_hf.M, lqrsol_hf.Qbar);
				var_J2 = St.LQRVar(x_mean, x_cov, U_lf(:, i), lqrsol_lf.Q, lqrsol_lf.S, lqrsol_lf.M, lqrsol_lf.Qbar);
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
				cost_hf_corr(i, :) = St.LQRCost(x0_rv_ext, Uopt_hf(:, i), lqrsol_hf);
				cost_lf_corr(i, :) = St.LQRCost(x0_rv_ext, Uopt_lf(:, i), lqrsol_lf);
			end
			
			corr = zeros(iters, iters);
			for i = 1:iters % HF iterations
				for j = 1:iters % LF iterations
					corr_mat = corrcoef(cost_hf_corr(i, :), cost_lf_corr(j, :));
					% we only need the cross correlation, diagnonal will be 1
					corr(i, j) = corr_mat(1,2); % rows are HF, cols are LF
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
				cost_hf(i,:) = St.LQRCost(x0_rv_ext, Uopt_hf(:, i), lqrsol_hf);
				cost_lf(i,:) = St.LQRCost(x0_rv_ext, Uopt_lf(:, i), lqrsol_lf);
			end
		end
		
	end
end