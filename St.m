% class that contains functions relevant to calculating statistics

classdef St
	methods(Static) % so that we don't need to create an object to use these functions
		
		% handles x0 having multiple samples with diag
		function cost = LQRCost(x0, lqrsol, U)
			Q = lqrsol.Q; S = lqrsol.S; M = lqrsol.M; Qbar = lqrsol.Qbar; Rbar = lqrsol.Rbar;
			cost = U'*(S'*Qbar*S + Rbar)*U + 2*x0'*M'*Qbar*S*U + diag(x0'*(M'*Qbar*M + Q)*x0);
		end
		
		% cost for multiple x0 samples with multiple Us
		function cost = LQRCostMulti(x0_rv, lqrsol, U) % TODO: U should be in rows
			items = size(U, 2);
			cost = zeros(items, size(x0_rv, 2)); % rows are items, cols are samples
			for i = 1:items
				cost(i, :) = St.LQRCost(x0_rv, lqrsol, U(:, i));
			end
		end
		
		function grad = LQRGrad(x0, lqrsol, U)
			S = lqrsol.S; M = lqrsol.M; Qbar = lqrsol.Qbar; Rbar = lqrsol.Rbar;
			% 2HU + 2q
			H = S'*Qbar*S + Rbar;
			q = x0'*M'*Qbar*S;
			grad = 2 * H * U + 2 * q';
		end
		
		function [f, g] = LQRCostwGrad(x0, lqrsol, U)
			f = St.LQRCost(x0, lqrsol, U);
			g = St.LQRGrad(x0, lqrsol, U);
		end
		
		function exp = LQRExp(x0_mean, x0_cov, lqrsol, U)
			Q = lqrsol.Q; S = lqrsol.S; M = lqrsol.M; Qbar = lqrsol.Qbar; Rbar = lqrsol.Rbar;
			K = U'*(S'*Qbar*S + Rbar)*U;
			L = 2 * M' * Qbar * S * U;
			N = M' * Qbar * M + Q;
			exp = K + x0_mean' * L + x0_mean' * N * x0_mean + trace(N * x0_cov);
		end
		
		function var = LQRVar(x0_mean, x0_cov, lqrsol, U)
			Q = lqrsol.Q; S = lqrsol.S; M = lqrsol.M; Qbar = lqrsol.Qbar;
			L = 2 * M' * Qbar * S * U;
			N = M' * Qbar * M + Q;
			var = L' * x0_cov * L + 2 * trace(N * x0_cov * N * x0_cov) + 4 * (x0_mean' * N + L') * x0_cov * N * x0_mean;
		end
		
		function cov = LQRCov(x0_mean, x0_cov, lqrsol_1, lqrsol_2, U_1, U_2)
			% K1 = Uopt_hf' * (lqrsol_hf.S' * lqrsol_hf.Qbar * lqrsol_hf.S + lqrsol_hf.Rbar) * Uopt_hf;
			L1 = 2 * lqrsol_1.M' * lqrsol_1.Qbar * lqrsol_1.S * U_1;
			N1 = lqrsol_1.M' * lqrsol_1.Qbar * lqrsol_1.M + lqrsol_1.Q;
			% K2 = Uopt_lf' * (lqrsol_lf.S' * lqrsol_lf.Qbar * lqrsol_lf.S + lqrsol_lf.Rbar) * Uopt_lf;
			L2 = 2 * lqrsol_2.M' * lqrsol_2.Qbar * lqrsol_2.S * U_2;
			N2 = lqrsol_2.M' * lqrsol_2.Qbar * lqrsol_2.M + lqrsol_2.Q;
			cov = L1'*x0_cov*L2 + 2*x0_mean'*N2*x0_cov*L1 + 2*x0_mean'*N1*x0_cov*L2 + 2*trace(N1*x0_cov*N2*x0_cov) + 4*x0_mean'*N1*x0_cov*N2*x0_mean;
		end
		
		function corr = LQRCorr(x0_mean, x0_cov, lqrsol_1, lqrsol_2, U_1, U_2)
			cov = St.LQRCov(x0_mean, x0_cov, lqrsol_1, lqrsol_2, U_1, U_2);
			var_J1 = St.LQRVar(x0_mean, x0_cov, lqrsol_1, U_1);
			var_J2 = St.LQRVar(x0_mean, x0_cov, lqrsol_2, U_2);
			corr = cov / sqrt(var_J1 * var_J2);
		end
		
		% correlation for multiple Us given the mean and cov of x0
		function corr = LQRCorrMulti(x0_mean, x0_cov, lqrsol_1, lqrsol_2, U_1, U_2)
			items = size(U_1, 2);
			corr = zeros(items, 1);
			for i = 1:items
				corr(i) = St.LQRCorr(x0_mean, x0_cov, lqrsol_1, lqrsol_2, U_1(:, i), U_2(:, i));
			end
		end
		
		% correlation for multiple Us across all pairs given the mean and cov of x0
		function corr = LQRCorrMulti2D(x0_mean, x0_cov, lqrsol_1, lqrsol_2, U_1, U_2)
			% TODO: Us should be in rows
			items = size(U_1, 2);
			corr = zeros(items, items);
			for i = 1:items
				for j = 1:items
					corr(i, j) = St.LQRCorr(x0_mean, x0_cov, lqrsol_1, lqrsol_2, U_1(:, i), U_2(:, j));
				end
			end
		end
		
		% correlation between the same rows of two matrices
		function corr = CorrMulti(cost_1, cost_2)
			items = size(cost_1, 1);
			corr = zeros(items, 1);
			for i = 1:items
				corr_mat = corrcoef(cost_1(i, :), cost_2(i, :));
				corr(i) = corr_mat(1, 2); % we only need the cross correlation, diagnonal will be 1
			end
		end
		
		% correlation between all pairs of rows of two matrices
		function corr = CorrMulti2D(cost_1, cost_2)
			items = size(cost_1, 1);
			corr = zeros(items, items);
			for i = 1:items % iterations in 1
				for j = 1:items % iterations in 2
					corr_mat = corrcoef(cost_1(i, :), cost_2(j, :));
					% rows are from 1, cols are from 2
					corr(i, j) = corr_mat(1, 2); % we only need the cross correlation, diagnonal will be 1
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