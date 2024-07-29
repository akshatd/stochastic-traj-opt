---
geometry: "margin=2cm"
---

# Multi-fidelity Stochastic Trajectory Optimization

## Overview

- Given the state space equations of a system, we can predict how it evolves over time with a given control input.
- Since we can predict the future states of the system, we can optimize the control inputs over time to make sure the system gets to a desired state (position and vecocity = 0 in this case).
- We also do not want to use too much control effort to get to the desired state.
- Putting both of these ideas into maths, we can construct a cost function that is dependent on how far away the state is to the desired state and how much control effort is being used.
- This function can then be minimized over the control inputs to see which control inputs can get the system to the desired state with minimal effort.
- Now some things in the system might not be known exactly, like the initial state of the system or the spring constant of a mass spring damper system.
- These unknowns can be modeled as random variables, and the cost function then becomes a random variable itself, since it is now a function of other random variables.
- This cost function can be estimated using Monte Carlo with random samples from the distribution of the parameters, which would give us an expected value of the cost function.
- If we use fewer samples in the Monte Carlo, the estimate will be less accurate and will be away from the true expected value of the cost function more often and vice-versa.
- This means that if we were to repeat the Monte Carlo simulation multiple times with a fixed number of samples, it will have some variance in the estimate of the cost function.
- Our goal then becomes to minimize this variance in the Monte Carlo estimate of the cost function given a finite total computational cost, so that every time we run Monte Carlo, we are sure that it wont be too far off from the true value.
- Now normally we would just increase the number of samples in the Monte Carlo to reduce the variance till we hit our computational budget.
- But what if we had access to a lower fidelity model of the system that could give us a rough estimate (maybe it has some bias in either the cost or the solution of control inputs) of the cost function with lower computational cost? This means we can run more samples in the lower fidelity model for the same computational cost!
- With this in mind, our goal of minimizing the variance in the Monte Carlo estimate of the cost function is now not just about maxxing out the number of samples, but finding a good balance between the the number of samples we run in the high and low fidelity models given a fixed total computational budget that minimizes the variance.
- If we have multiple fidelity models, we will have to determine how to distribute the samples among the models.
- So to summarize
  - We have a system, whose trajectory we want to optimize using a cost function
  - Some parameters in the system are random, so the cost function is also random
  - To get an expectation of the cost function, we use Monte Carlo using random samples of the random parameters
  - If we have a low fidelity model as well, we can run a lot of samples in that model or run a few samples in the high fidelity model.
  - We want the Monte Carlo estimate to be close to the true value of the cost function, so we want to minimize its run-to-run variance.
  - Hence we try to distribute the samples among the models for a single run of the Monte Carlo so its variance is minimized.

## Dynamics

We have a simple mass spring damper system:

- Mass: m (chosen as 5)
- Spring constant: k (chosen as 2)
- Damping coefficient: b (chosen as 0.5)
- External force: F
- Displacement: x
- Velocity: v = $\dot{x}$
- Acceleration: a = $\ddot{x}$
- Force on Mass = $ma = m\ddot{x}$

Hence we have the equality

$$
\begin{aligned}
m\ddot{x} &= F - b\dot{x} - kx \\
\implies \ddot{x} &= \frac{F}{m} - \frac{b}{m}\dot{x} - \frac{k}{m}x
\end{aligned}
$$

Let the state of the dynamics be $[x, \dot{x}]$ where the control input $u=F$. We just observe the position $x$. Hence we have the state space equations

$$
\begin{aligned}
\begin{bmatrix}
\dot{x} \\
\ddot{x}
\end{bmatrix} &= \begin{bmatrix}
0 & 1 \\
-\frac{k}{m} & -\frac{b}{m}
\end{bmatrix} \begin{bmatrix}
x \\
\dot{x}
\end{bmatrix} + \begin{bmatrix}
0 \\
\frac{1}{m}
\end{bmatrix} u \\
y &=
\begin{bmatrix}
1 & 0
\end{bmatrix} \begin{bmatrix}
x \\
\dot{x}
\end{bmatrix}
\end{aligned}
$$

Thus, we have a linear system in the form $\dot{x} = Ax + Bu$ and $y = Cx$. Where

$$
\begin{aligned}
A &= \begin{bmatrix}
0 & 1 \\
-\frac{k}{m} & -\frac{b}{m}
\end{bmatrix} \\
B &= \begin{bmatrix}
0 \\
\frac{1}{m}
\end{bmatrix} \\
C &= \begin{bmatrix}
1 & 0
\end{bmatrix}
\end{aligned}
$$

This linear system will be converted to discrete time using MATLAB to give us the form $x_{k+1} = A x_k + B u_k$ and $y_k = C x_k$ for a chosen sampling time $T_s$.

![Open Loop Dynamics of the Mass Spring Damper System](figs/ol_det.svg)

## Deterministic LQR

To solve the LQR problem with reference tracking, we need to define the cost function to be minimzed. This is defined as

$$
\min_{u_0, u_1 ... u_{N-1}}J = \sum_{k=0}^{N-1} (e_k^T Q_e e_k + \Delta u_k^T R \Delta u_k)
$$

Where:

- $u_, u_1 ... u_{N-1}$ are the control inputs at each timestep of the trajectory
- $J$ is the total cost to be minimized
- $N$ is the time horizon
- $x_k$ is the state vector at time step $k$
- $u_k$ is the control input vector at time step $k$
- $\Delta u_k = u_k - u_{k-1}$ is the change in control input at time step $k$
- $r_k$ is the reference output at time step $k$
- $e_k = y_k - r_k = Cx - r_k$ is the error in the output at time step $k$
- $Q_e$ is the error cost matrix
- $R$ is the control cost matrix
<!-- - $P$ is the terminal cost matrix -->

Here we extend the usual state with other variables to make it easier to track the error in the output and the change in control input. The state vector is now

$$
x_k^{ext} = \begin{bmatrix}
x_k \\
u_{k-1} \\
r_k
\end{bmatrix}
$$

And the state space equations are

$$
\begin{aligned}
x_{k+1}^{ext} &= \begin{bmatrix}
A & B & 0 \\
0 & \mathbb{I}_{n_u \times n_u} & 0 \\
0 & 0 & \mathbb{I}_{n_r \times n_r}
\end{bmatrix} \begin{bmatrix}
x_k \\
u_{k-1} \\
r_k
\end{bmatrix} + \begin{bmatrix}
B \\
\mathbb{I}_{n_u \times n_u} \\
0
\end{bmatrix} \Delta u_k \\
\end{aligned}
$$

with

$$
\begin{aligned}
A^{ext} &= \begin{bmatrix}
A & B & 0 \\
0 & \mathbb{I}_{n_u \times n_u} & 0 \\
0 & 0 & \mathbb{I}_{n_r \times n_r}
\end{bmatrix} \\
B^{ext} &= \begin{bmatrix}
B \\
\mathbb{I}_{n_u \times n_u} \\
0
\end{bmatrix}
\end{aligned}
$$

Stacking the state and control input vectors, we have

$$
\begin{aligned}
X &= \begin{bmatrix}
x_1^{ext} \\
x_2^{ext} \\
\vdots \\
x_N^{ext}
\end{bmatrix} \text{with size } N n_x \times 1 \\
U &= \begin{bmatrix}
\Delta u_0 \\
\Delta u_1 \\
\vdots \\
\Delta u_{N-1}
\end{bmatrix} \text{with size } N n_u \times 1
\end{aligned}
$$

Additionally, The matrix $Q$ for the extended state should only penalize the error, so for the extended state, the Q matrix can be derived as follows

$$
\begin{aligned}
e_k^T Q_e e_k &= (y_k - r_k)^T Q_e (y_k - r_k) \\
&= \left(\begin{bmatrix} C & 0 & -\mathbb{I_{n_r \times n_r}} \end{bmatrix}
\begin{bmatrix}
x_k \\
u_{k-1} \\
r_k
\end{bmatrix}\right)^T Q_e
\left(\begin{bmatrix} C & 0 & -\mathbb{I_{n_r \times n_r}} \end{bmatrix}
\begin{bmatrix}
x_k \\
u_{k-1} \\
r_k
\end{bmatrix}\right) \\
&= x_k^{ext^T} \begin{bmatrix} C & 0 & -\mathbb{I_{n_r \times n_r}} \end{bmatrix}^T Q_e
\begin{bmatrix} C & 0 & -\mathbb{I_{n_r \times n_r}} \end{bmatrix} x_k^{ext} \\
\implies Q &= \begin{bmatrix} C & 0 & -\mathbb{I_{n_r \times n_r}} \end{bmatrix}^T Q_e
\begin{bmatrix} C & 0 & -\mathbb{I_{n_r \times n_r}} \end{bmatrix}
\end{aligned}
$$

I will represent the new extended state $x_k^{ext}$ and state transition matrices $A^{ext}$, $B^{ext}$ as just $x_k$, $A$, $B$ to make things less cluttered. State transition formula for a linear system is $x_{k+1} = A x_k + B u_k$. Hence, starting from $x_0$, we can write the state vector at any time step as

$$
\begin{aligned}
x_1 &= A x_0 + B u_0 \\
x_2 &= A x_1 + B u_1 = A(A x_0 + B u_0) + B u_1 = A^2 x_0 + A B u_0 + B u_1 \\
\vdots \\
x_N &= A x_{N-1} + B u_{N-1} = A^{N-1} x_0 + A^{N-2} B u_0 + A^{N-3} B u_1 + \ldots + B u_{N-1} \\
\implies x_k &= A^k x_0 + \sum_{i=0}^{k-1} A^{k-1-i} B u_i
\end{aligned}
$$

Using this formula, we can construct matrices $S$ and $M$ such that $X = S U + M x_0$. The matrices $S$ and $M$ are given by

$$
\begin{aligned}
S &= \begin{bmatrix}
B & 0 & \ldots & 0 \\
AB & B & \ldots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
A^{N-1} B & A^{N-2} B & \ldots & B
\end{bmatrix} \text{with size } Nn_x \times Nn_u \\
M &= \begin{bmatrix}
A \\
A^2 \\
\vdots \\
A^N
\end{bmatrix} \text{with size } Nn_x \times n_x
\end{aligned}
$$

The cost matrices can be stacked in a similar way to get $\bar{Q}$ and $\bar{R}$

$$
\begin{aligned}
\bar{Q} &= \begin{bmatrix}
Q & 0 & \ldots & 0 \\
0 & Q & \ldots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \ldots & 0
\end{bmatrix} \text{with size } Nn_x \times Nn_x \\
\bar{R} &= \begin{bmatrix}
R & 0 & \ldots & 0 \\
0 & R & \ldots & 0 \\
\vdots & \vdots & \ddots & \vdots \\
0 & 0 & \ldots & R
\end{bmatrix} \text{with size } Nn_u \times Nn_u
\end{aligned}
$$

Note that the final state is not being penalized with the terminal cost matrix $P$ as we do not want to penalize the state for being far away from the origin.

The cost function in terms of the new matrices $X$, $U$, $S$, $M$, $\bar{Q}$ and $\bar{R}$ is

$$
\begin{aligned}
J &= x_N^T P x_N + \sum_{k=0}^{N-1} (x_k^T Q x_k + u_k^T R u_k) \\
&= X^T \bar{Q} X + U^T \bar{R} U + x_0^T Q x_0 \text{ (extra $x_0$ term because it is not in $X$)} \\
&= (S U + M x_0)^T \bar{Q} (S U + M x_0) + U^T \bar{R} U + x_0^T Q x_0 \\
% & \text{note } (A + B)^T = A^T + B^T \text{ and } (AB)^T = B^T A^T \\
&= (U^T S^T + x_0^T M^T) \bar{Q} (S U + M x_0) + U^T \bar{R} U + x_0^T Q x_0 \\
% & \text{note } A(B + C) = AB + AC \\
&= (U^T S^T + x_0^T M^T) (\bar{Q} S U + \bar{Q} M x_0) + U^T \bar{R} U + x_0^T Q x_0 \\
&= U^T S^T \bar{Q} S U + U^T S^T \bar{Q} M x_0 + x_0^T M^T \bar{Q} S U + x_0^T M^T \bar{Q} M x_0 + U^T \bar{R} U + x_0^T Q x_0 \\
% & \text{note } \bar{Q} \text{ is symmetric so it doesnt need to be transposed} \\
&= U^T S^T \bar{Q} S U + (x_0^T M^T \bar{Q} S U)^T + x_0^T M^T \bar{Q} S U + x_0^T M^T \bar{Q} M x_0 + U^T \bar{R} U + x_0^T Q x_0 \\
% & \text{note } x_0^T M^T \bar{Q} S U \text{ is a scalar(check its dimensions), so its transpose is itself} \\
&= U^T S^T \bar{Q} S U + x_0^T M^T \bar{Q} S U + x_0^T M^T \bar{Q} S U + x_0^T M^T \bar{Q} M x_0 + U^T \bar{R} U + x_0^T Q x_0 \\
% & \text{combining all quadratic and repeated terms} \\
&= U^T (S^T \bar{Q} S + \bar{R}) U + 2 x_0^T M^T \bar{Q} S U + x_0^T(M^T \bar{Q} M + Q) x_0\\
\end{aligned}
$$

Let

$$
\begin{aligned}
H &= S^T \bar{Q} S + \bar{R} = H^T \\
&\text{(quadratic multiplication by diagnonal($\bar{Q}$) results in a symmetric, R is symmetric)} \\
q &= (x_0^T M^T \bar{Q} S)^T = S^T \bar{Q} M x_0 \\
c &= x_0^T (M^T \bar{Q} M + Q) x_0
\end{aligned}
$$

The cost function can now be written as a quadratic function of the control inputs $U$

$$
J = U^T H U + 2 q^T U + c
$$

Which has a gradient with respect to $U$ as

$$
\Delta J = 2 H U + 2 q
$$

Assuming minimum exists, we can set the gradient to zero to get the optimal control input $U^*$

$$
\begin{aligned}
\Delta J &= 2 H U + 2 q = 0 \\
\implies U^* &= -H^{-1} q
\end{aligned}
$$

![Closed Loop with LQR Control](figs/0.05_cl_det.svg){width=50%}
![Closed Loop with LQR Control](figs/0.5_cl_det.svg){width=50%}

## Stochastic LQR

We will consider the case where the initial state $x_0$ is stochastic, with mean $\mu$ and covariance $\Sigma$.
The control input vector $U$ is still deterministic.

In this case, the state vector $X$ is stochastic, owing to the the initial condition $x_0$ being a vector of random variables. The state propagation equations are now expectations

$$
\begin{aligned}
\mathbb{E}[x_1] &= \mathbb{E}[A x_0 + B u_0] = A \mathbb{E}[x_0] + B u_0 \\
\mathbb{E}[x_2] &= \mathbb{E}[A x_1 + B u_1] = A \mathbb{E}[x_1] + B u_1 = A^2 \mathbb{E}[x_0] + A B u_0 + B u_1 \\
\vdots \\
\mathbb{E}[x_N] &= \mathbb{E}[A x_{N-1} + B u_{N-1}] = A \mathbb{E}[x_{N-1}] + B u_{N-1} \\
 &= A^N \mathbb{E}[x_0] + A^{N-1} B u_0 + A^{N-2} B u_1 + \ldots + B u_{N-1} \\
\end{aligned}
$$

![Open Loop Dynamics with Stochastic Initial Conditions](figs/ol_stoch_init.svg)

The matrices $S$ and $M$ remain the same, and we have the following expression for the expectation of the state vector $X$

$$
\mathbb{E}[X] = S U + M \mathbb{E}[x_0]
$$

The cost matrices $\bar{Q}$ and $\bar{R}$ are also the same. The cost function is now a random variable $\hat{J}$, with an expectation and variance.

$$
\hat{J} \sim \mathcal{N}(\mathbb{E}[\hat{J}], Var[\hat{J}])
$$

The expected value of the cost function over all possible realizations of the initial conditions is

$$
\begin{aligned}
\mathbb{E}[\hat{J}] &= \mathbb{E}[x_N^T P x_N + \sum_{k=0}^{N-1} (x_k^T Q x_k + u_k^T R u_k)] \\
&= \mathbb{E}[X^T \bar{Q} X + U^T \bar{R} U + x_0^T Q x_0] \\
&= \mathbb{E}[X^T \bar{Q} X] + \mathbb{E}[U^T \bar{R} U] + \mathbb{E}[x_0^T Q x_0] \\
&= \mathbb{E}[(S U + M x_0)^T \bar{Q} (S U + M x_0)] + U^T \bar{R} U + \mathbb{E}[x_0^T Q x_0] \\
% & \text{quadratic form: https://en.wikipedia.org/wiki/Quadratic_form_(statistics)} \\
&= \mathbb{E}[U^T S^T \bar{Q} S U + U^T S^T \bar{Q} M x_0 + x_0^T M^T \bar{Q} S U + x_0^T M^T \bar{Q} M x_0] + U^T \bar{R} U + \mathbb{E}[x_0]^T Q \mathbb{E}[x_0] + \text{tr}(Q \Sigma) \\
&= U^T S^T \bar{Q} S U + \mathbb{E}[U^T S^T \bar{Q} M x_0] + \mathbb{E}[x_0^T M^T \bar{Q} S U] + \mathbb{E}[x_0^T M^T \bar{Q} M x_0] + U^T \bar{R} U + \mathbb{E}[x_0]^T Q \mathbb{E}[x_0] + \text{tr}(Q \Sigma) \\
&= U^T (S^T \bar{Q} S + \bar{R}) U + U^T S^T \bar{Q} M \mathbb{E}[x_0] + \mathbb{E}[x_0^T] M^T \bar{Q} S U + \mathbb{E}[x_0^T M^T \bar{Q} M x_0] + \mathbb{E}[x_0]^T Q \mathbb{E}[x_0] + \text{tr}(Q \Sigma) \\
&= U^T (S^T \bar{Q} S + \bar{R}) U + \mathbb{E}[x_0^T] M^T \bar{Q} S U + \mathbb{E}[x_0^T] M^T \bar{Q} S U + \mathbb{E}[x_0]^T M^T \bar{Q} M \mathbb{E}[x_0] + \text{tr}(M^T \bar{Q} M \Sigma) + \mathbb{E}[x_0]^T Q \mathbb{E}[x_0] + \text{tr}(Q \Sigma) \\
&= U^T (S^T \bar{Q} S + \bar{R}) U + 2 \mathbb{E}[x_0^T] M^T \bar{Q} S U + \mathbb{E}[x_0]^T (M^T \bar{Q} M + Q) \mathbb{E}[x_0] + \text{tr}(M^T \bar{Q} M \Sigma)  + \text{tr}(Q \Sigma) \\
\end{aligned}
$$

Now assuming $\mathbb{E}[x_0]$ is calculated using a Monte Carlo estimator, we can write the cost function as

$$
\begin{aligned}
H &= S^T \bar{Q} S + \bar{R} = H^T \\
q &= (\mathbb{E}[x_0]^T M^T \bar{Q} S)^T = S^T \bar{Q} M \mathbb{E}[x_0] \\
c &= \mathbb{E}[x_0]^T (M^T \bar{Q} M + Q) \mathbb{E}[x_0] + \text{tr}(M^T \bar{Q} M \Sigma)  + \text{tr}(Q \Sigma)
\end{aligned}
$$

The cost function can now be written as a quadratic function of the control inputs $U$

$$
J = U^T H U + 2 q^T U + c
$$

Which has a gradient with respect to $U$ as

$$
\Delta J = 2 H U + 2 q
$$

Assuming minimum exists, we can set the gradient to zero to get the optimal control input $U^*$

$$
\begin{aligned}
\Delta J &= 2 H U + 2 q = 0 \\
\implies U^* &= -H^{-1} q
\end{aligned}
$$

Which proves that the optimal control input when the initial state is stochastic is the same as when the initial state is deterministic if that the mean of the stochastic initial state is the same as the deterministic initial state given enough samples.

To visualize this, we can optmize individual trajectories with randomly sampled initial conditions and plot the average trajectory. The average trajectory is the solution to the stochastic LQR problem.

![Closed Loop with Stochastic LQR Control](figs/0.05_cl_stoch_init.svg)
![Closed Loop with Stochastic LQR Control](figs/0.5_cl_stoch_init.svg)

To calculate the variance, we can rewrite the expression of the cost function

$$
\begin{aligned}
J &= U^T(S^T \bar{Q} S + \bar{R}) U + 2 x_0^T M^T \bar{Q} S U + x_0^T(M^T \bar{Q} M + Q) x_0 \\
\text{Let}& \\
K &= U^T(S^T \bar{Q} S + \bar{R}) U \\
L &= 2 M^T \bar{Q} S U \\
N &= M^T \bar{Q} M + Q \\
\text{Then}& \\
J &= K + x_0^T L + x_0^T N x_0 \\
\end{aligned}
$$

Similarly, expectation of the random cost function can be rewritten as

$$
\begin{aligned}
\mathbb{E}[\hat{J}] &= U^T(S^T \bar{Q} S + \bar{R}) U + 2 \mathbb{E}[x_0^T] M^T \bar{Q} S U + \mathbb{E}[x_0]^T(M^T \bar{Q} M + Q) \mathbb{E}[x_0] + \text{tr}(M^T \bar{Q} M \Sigma) + \text{tr}(Q \Sigma) \\
\text{Let}& \\
O &= \text{tr}(M^T \bar{Q} M \Sigma) + \text{tr}(Q \Sigma) \\
\text{Then}& \\
\mathbb{E}[\hat{J}] &= K + \mathbb{E}[x_0^T] L + \mathbb{E}[x_0]^T N \mathbb{E}[x_0] + O \\
\end{aligned}
$$

Then, the variance can be written as

$$
\begin{aligned}
Var[\hat{J}] &= \mathbb{E}[(\hat{J} - \mathbb{E}[\hat{J}])^2] \\
&= \mathbb{E}[(K + x_0^T L + x_0^T N x_0 - (K + \mathbb{E}[x_0^T] L + \mathbb{E}[x_0]^T N \mathbb{E}[x_0] + O))^2] \\
&= \mathbb{E}[(x_0^T L + x_0^T N x_0 - \mathbb{E}[x_0^T] L - \mathbb{E}[x_0]^T N \mathbb{E}[x_0] - O)^2] \\
&= \mathbb{E}[(x_0^T L - \mathbb{E}[x_0^T] L + x_0^T N x_0 - \mathbb{E}[x_0]^T N \mathbb{E}[x_0] - O)^2] \\
&= \mathbb{E}[((x_0^T - \mathbb{E}[x_0]^T) L + (x_0 + \mathbb{E}[x_0])^T N (x_0 - \mathbb{E}[x_0]) + O)^2] \\
\end{aligned}
$$

Here we use a simple Monte Carlo estimator to estimate the expected value of the cost function over some realizations of the initial conditions. The variance of this cost function is varies with the number of samples used in the Monte Carlo estimator.

![Variance of the Cost Function with Number of Samples](figs/0.05_mc_variance.svg){width=50%}
![Variance of the Cost Function with Number of Samples](figs/0.5_mc_variance.svg){width=50%}
