import numpy as np
import scipy as sp


class LQRMSD:
    def __init__(self, t_step, t_horiz, Qe, Re, Pe):

        # System params
        m = 5  # mass
        k = 2  # spring coefficient
        c = 0.5  # damping coefficient
        Ac = np.array([
            [0, 1],
            [-k / m, -c / m]
        ])
        Bc = np.array([
            [0],
            [1 / m]
        ])
        Cc = np.array([
            [1, 0],
            [0, 1]
        ])
        Dc = np.array([
            [0],
            [0]
        ])

        self.t_step = t_step

        # Discretize the system
        sys_cont = sp.signal.lti(Ac, Bc, Cc, Dc)
        sys_disc = sys_cont.to_discrete(t_step, method='zoh')
        self.A, self.B, self.Q, self.R, self.P = self.extendState(sys_disc, Qe, Re, Pe)
        self.nx = self.A.shape[0]  # Number of states
        self.nu = self.B.shape[1]  # Number of inputs

        # Set up LQR problem
        self.N = int(t_horiz / t_step)  # Prediction horizon
        self.Kopt, self.S, self.M, self.Qbar, self.Rbar = self.solveLQR()

    def obj(self, x0, U):
        """
        Compute the cost for the LQR problem given the multiple x0 and the control input U.
        """
        x0 = x0.reshape(self.nx, -1)  # fix shape
        U = U.reshape(self.N*self.nu, -1)  # fix shape
        Q, S, M, Qbar, Rbar = self.Q, self.S, self.M, self.Qbar, self.Rbar
        # Compute the cost
        cost = U.T @ (S.T @ Qbar @ S + Rbar) @ U + 2 * x0.T @ M.T @ Qbar @ S @ U + \
            np.diag(x0.T @ (M.T @ Qbar @ M + Q) @ x0).reshape(-1, 1)
        return cost

    def grad(self, x0, U):
        """
        Compute the gradient of the cost for the LQR problem given the multiple x0 and the control input U.
        """
        x0 = x0.reshape(self.nx, -1)  # fix shape
        U = U.reshape(self.N*self.nu, -1)  # fix shape
        S, M, Qbar, Rbar = self.S, self.M, self.Qbar, self.Rbar
        # Compute the gradient
        H = S.T @ Qbar @ S + Rbar
        q = x0.T @ M.T @ Qbar @ S
        grad = 2 * H @ U + 2 * q.T
        return grad

    def exp(self, x0_mean, x0_cov, U):
        """
        Compute the expected cost for the LQR problem given the mean and covariance of x0 and the control input U.
        """
        Q, S, M, Qbar, Rbar = self.Q, self.S, self.M, self.Qbar, self.Rbar
        # Compute the expected cost
        K = U.T @ (S.T @ Qbar @ S + Rbar) @ U
        L = 2 * M.T @ Qbar @ S @ U
        N = M.T @ Qbar @ M + Q
        exp = K + x0_mean.T @ L + x0_mean.T @ N @ x0_mean + np.trace(N @ x0_cov)
        return exp

    def var(self, x0_mean, x0_cov, U):
        """
        Compute the variance of the LQR problem given the mean and covariance of x0 and the control input U.
        """
        Q, S, M, Qbar = self.Q, self.S, self.M, self.Qbar
        # Compute the variance
        L = 2 * M.T @ Qbar @ S @ U
        N = M.T @ Qbar @ M + Q
        var = L.T @ x0_cov @ L + 2 * np.trace(N @ x0_cov @ N @ x0_cov) + \
            4 * (x0_mean.T @ N + L.T) @ x0_cov @ N @ x0_mean
        return var

    def extendState(self, system, Qe, Re, Pe):
        A = system.A
        B = system.B
        C = system.C

        nx = A.shape[0]  # Number of states
        nu = B.shape[1]  # Number of inputs
        ny = C.shape[0]  # Number of measurements
        nr = ny  # Number of references (equivalent to measurments)

        # Extend states are [x_k, u_k-1, r_k]
        Aext = np.block([
            [A, B, np.zeros((nx, nr))],
            [np.zeros((nu, nx)), np.eye(nu), np.zeros((nu, nr))],
            [np.zeros((nr, nx)), np.zeros((nr, nu)), np.eye(nr)]
        ])

        Bext = np.vstack([B, np.eye(nu), np.zeros((nr, nu))])

        E = np.hstack([C, np.zeros((ny, nu)), -np.eye(nr)])

        Q = E.T @ Qe @ E
        R = Re
        P = sp.linalg.block_diag(Pe, np.zeros((Aext.shape[0] - nx, Aext.shape[0] - nx)))

        return Aext, Bext, Q, R, P

    def solveLQR(self):
        # Compute the matrices S, M, Qbar, Rbar, and K0N
        # for the unconstrained LQ-MPC problem
        #
        # Inputs:
        #   N: Prediction horizon
        #   A: State transition matrix
        #   B: Input matrix
        #   Q: State cost matrix
        #   R: Input cost matrix
        #   P: Terminal state cost matrix

        nx = self.nx
        nu = self.nu
        N = self.N
        A = self.A
        B = self.B
        Q = self.Q
        R = self.R
        P = self.P

        # Initialize matrices
        S = np.zeros((N * nx, N * nu))
        M = np.zeros((N * nx, nx))
        Qbar = np.zeros((N * nx, N * nx))
        Rbar = np.zeros((N * nu, N * nu))

        # Compute the first column of S
        for i in range(1, N + 1):
            rowStart = (i - 1) * nx
            rowEnd = i * nx
            S[rowStart:rowEnd, :nu] = sp.sparse.linalg.matrix_power(A, i - 1) @ B

        # Pad the first column and set it to other columns of S
        for i in range(2, N + 1):
            colStart = (i - 1) * nu
            colEnd = i * nu
            zeroRows = (i - 1) * nx
            zeroCols = nu
            S[:, colStart:colEnd] = np.vstack(
                [np.zeros((zeroRows, zeroCols)), S[0:-zeroRows, :nu]])

        # Compute first row of M
        M[0:nx, :] = A
        # Compute the rest of M
        for i in range(2, N + 1):
            rowStart = (i - 1) * nx
            rowEnd = i * nx
            # just multiply the previous rows by A to get higher powers
            M[rowStart:rowEnd, :] = A @ M[rowStart - nx:rowEnd - nx, :]

        # Compute Qbar except for the last row
        for i in range(1, N + 1):
            # Q is square so we can reuse indices
            rowStart = (i - 1) * nx
            rowEnd = i * nx
            temp = Q
            if i == N:
                temp = P
            Qbar[rowStart:rowEnd, rowStart:rowEnd] = temp

        # Compute Rbar
        for i in range(1, N + 1):
            # R is square so we can reuse indices
            rowStart = (i - 1) * nu
            rowEnd = i * nu
            Rbar[rowStart:rowEnd, rowStart:rowEnd] = R
        # Compute Optimal Control Gain
        Kopt = -np.linalg.inv(S.T @ Qbar @ S + Rbar) @ S.T @ Qbar @ M
        return Kopt, S, M, Qbar, Rbar
