## Dynamics

The objective is to solve a reference tracking problem for the position and velocity of a simple mass-spring-damper system.
The mass-spring-damper system has the following parameters:

- $m = 5$ g : Mass
- $k = 2$ N/m : Spring constant
- $c = 0.5$ Ns/m : Damping coefficient

Given an external force $F$ and the acceleration of the system $\ddot{x}$,
the dynamics of the system can be written as the equality of forces acting on the mass:

$$
\begin{aligned}
   m\ddot{x} &= F - c\dot{x} - kx \\
   \implies \ddot{x} &= \frac{F}{m} - \frac{c}{m}\dot{x} - \frac{k}{m}x
\end{aligned}
$$

We have the state as $[x \dot{x}]^T$ containing the position and velocity of the mass, both of which are observed.
The control input $u=F$ is applied to the mass, yielding the following state space model:

$$
\begin{aligned}
  \begin{bmatrix}
    \dot{x} \\
    \ddot{x}
  \end{bmatrix} &=
  \begin{bmatrix}
    0 & 1 \\
    -\frac{k}{m} &  -\frac{c}{m}
  \end{bmatrix}
  \begin{bmatrix}
    x \\
    \dot{x}
  \end{bmatrix} +
  \begin{bmatrix}
    0 \\
    \frac{1}{m}
  \end{bmatrix} u \\
  y &=
  \begin{bmatrix}
    1 & 0 \\
    0 & 1
  \end{bmatrix}
  \begin{bmatrix}
    x \\
    \dot{x}
  \end{bmatrix}
\end{aligned}
$$

Thus, we have a standard continuous-time LTI system $\dot{x} = A_cx + B_cu$ and $y = C_cx$, where

$$
\begin{aligned}
   A_c &=
   \begin{bmatrix}
   0 & 1 \\
   -\frac{k}{m} & -\frac{c}{m}
   \end{bmatrix} \\
   B_c &=
   \begin{bmatrix}
   0 \\
   \frac{1}{m}
   \end{bmatrix} \\
   C_c &=
   \begin{bmatrix}
   1 & 0 \\
   0 & 1
   \end{bmatrix}
\end{aligned}
$$

## LQR

The continuous-time system can be converted to discrete time using such that $x_{k+1} = A x_k + B u_k$ and $y_k = C x_k$ for a chosen sampling time $T_s$.
The $A$ and $B$ matrices for the discrete time system are then used to formulate the cost function for the LQR problem.

We set up an LQR controller to track the reference position and velocity of the mass-spring-damper system.

The LQR controller minimizes the cost function of the extended state containing the reference output $r_k$

$$
x_k^{ext} =
\begin{bmatrix}
  x_k \\
  u_{k-1} \\
  r_k
\end{bmatrix}
$$

With the state space model

$$
\begin{aligned}
  x_{k+1}^{ext} &=
  \begin{bmatrix}
    A & B & 0 \\
    0 & \mathbb{I}_{n_u \times n_u} & 0 \\
    0 & 0 & \mathbb{I}_{n_r \times n_r}
  \end{bmatrix}
  \begin{bmatrix}
    x_k \\
    u_{k-1} \\
    r_k
  \end{bmatrix} +
  \begin{bmatrix}
    B \\
    \mathbb{I}_{n_u \times n_u} \\
    0
  \end{bmatrix} \Delta u_k \\
\end{aligned}
$$

And the cost function

$$
J = \sum_{k=0}^{N-1} x_k^{extT} Q x^{ext}_k + \Delta u_k^T R \Delta u_k
$$

We construct the matrices $X$, $U$, $S$ and $M$ such that $X = SU + Mx_0$.

$$
\begin{aligned}
  X &=
  \begin{bmatrix}
    x_1^{ext} \\
    x_2^{ext} \\
    \vdots \\
    x_N^{ext}
  \end{bmatrix}
  U &=
  \begin{bmatrix}
    \Delta u_0 \\
    \Delta u_1 \\
    \vdots \\
    \Delta u_{N-1}
  \end{bmatrix}
  S &=
  \begin{bmatrix}
    B & 0 & \ldots & 0 \\
    AB & B & \ldots & 0 \\
    \vdots & \vdots & \ddots & \vdots \\
    A^{N-1} B & A^{N-2} B & \ldots & B
  \end{bmatrix}
  M &=
  \begin{bmatrix}
    A \\
    A^2 \\
    \vdots \\
    A^N
  \end{bmatrix}
\end{aligned}
$$

The cost function in terms of the new matrices $X$, $U$, $S$, $M$, $\bar{Q}$ and $\bar{R}$ is

$$
\begin{aligned}
   J &= \sum_{k=0}^{N-1} x_k^T Q x_k + \Delta u_k^T R \Delta u_k \\
   % \text{ (extra $x_0$ term because it is not in $X$)} \\
   &= X^T \bar{Q} X + U^T \bar{R} U + x_0^T Q x_0 \\
   &= (S U + M x_0)^T \bar{Q} (S U + M x_0) + U^T \bar{R} U + x_0^T Q x_0 \\
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
\frac{\partial J}{\partial U} = 2 H U + 2 q^T
$$
