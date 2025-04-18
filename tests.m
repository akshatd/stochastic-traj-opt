clc; close all; clear;

% setup
nx = 3;
mu = 1.21;
sigma = 1;
x_mean = repelem(mu, nx)';
x_cov = sigma * eye(nx);
A = randn(nx, 1);
B = randn(nx, nx);
B = B * B'; % make it symmetric
samples = 10000;
x = mvnrnd(x_mean, x_cov, samples)'; % transpose to make it nx by samples
% change the vars to represent reality
x_mean = mean(x, 2);
x_cov = cov(x');

% expected values
disp("Expected values");

% test constant
disp("Constant");
st_mean = mean(x, 2);
an_mean = x_mean;
fprintf("Analytic: %f, Statistic: %f\n", an_mean, st_mean);

% test linear
disp("Linear");
term = x' * A;
st_mean = mean(term);
an_mean = x_mean' * A;
fprintf("Analytic: %f, Statistic: %f\n", an_mean, st_mean);

% test quadratic
disp("Quadratic");
term = diag(x' * B * x); % diag cos it contains the square term of x with itself
st_mean = mean(term);
an_mean = trace(B * x_cov) + x_mean' * B * x_mean;
fprintf("Analytic: %f, Statistic: %f\n", an_mean, st_mean);

% test cubic
disp("Cubic");
term = 0;
for i = 1:samples
	term = term + (x(:, i)' * A) * (x(:, i)' * B * x(:, i));
end
st_mean = term/samples;
% chatgpt, doesnt work
% an_mean = (x_mean' * B * x_mean) * (x_mean' * A) + trace(B * x_cov) * x_mean' * A
% UMGPT, works
% an_mean = 2 * x_mean' * B * x_cov * A + trace(B * x_cov) * x_mean' * A + (x_mean' * B * x_mean) * (x_mean' * A);
an_mean = (2 * x_mean' * B * x_cov + trace(B * x_cov) * x_mean' + x_mean' * B * x_mean * x_mean') * A;
fprintf("Analytic: %f, Statistic: %f\n", an_mean, st_mean);

% test quartic
disp("Quartic");
term = 0;
for i = 1:samples
	term = term + (x(:, i)' * B * x(:, i)) * (x(:, i)' * B * x(:, i));
end
st_mean = term/samples;
% textbook, no simplification, works
% an_mean = trace(B * x_cov* (B + B') * x_cov) + x_mean' * (B + B') * x_cov * (B + B') * x_mean + (trace(B * x_cov) + x_mean' * B * x_mean) * (trace(B * x_cov) + x_mean' * B * x_mean)
% mine, simplified, works
an_mean = 2 * trace(B * x_cov * B * x_cov) + 4 * x_mean' * B * x_cov * B * x_mean + (trace(B * x_cov) + x_mean' * B * x_mean)^2;
% chatgpt, doesnt work
% an_mean = (x_mean' * B * x_mean)^2 + 4 * (x_mean' * B * x_mean) * trace(B * x_cov) + 2 * trace(B * x_cov * B * x_cov) + trace(B * x_cov)^2
% UMGPT, doesnt work
% an_mean = (x_mean' * B * x_mean)^2 + 2 * (x_mean' * B * x_mean) * trace(B * x_cov) + trace(B * x_cov)^2 + 2 * trace(B * x_cov * B * x_cov)
fprintf("Analytic: %f, Statistic: %f\n", an_mean, st_mean);

% variances
disp(newline + "Variances");

% test linear
disp("Linear");
term = x' * A;
st_var = var(term);
an_var = A' * x_cov * A;
fprintf("Analytic: %f, Statistic: %f\n", an_var, st_var);

% test quadratic
disp("Quadratic");
term = diag(x' * B * x);
st_var = var(term);
an_var = trace(B * x_cov * (B + B') * x_cov) + x_mean' * (B + B') * x_cov * (B + B') * x_mean;
fprintf("Analytic: %f, Statistic: %f\n", an_var, st_var);

% test combined variance
disp("Combined");
term = x' * A + diag(x' * B * x);
st_var = var(term);
an_var = 2 * trace(B * x_cov * B * x_cov) + A' * x_cov * A + 4 * (x_mean' * B + A') * x_cov * B * x_mean;
fprintf("Analytic: %f, Statistic: %f\n", an_var, st_var);

% covariance
disp(newline + "Covariance");
K1 = randn;
K2 = randn;
L1 = randn(nx, 1);
L2 = randn(nx, 1);
N1 = randn(nx, nx);
N1 = N1 * N1';
N2 = randn(nx, nx);
N2 = N2 * N2';
term1 = zeros(1, samples);
term2 = zeros(1, samples);
for i = 1:samples
	term1(i) = K1 + x(:, i)' * L1 + x(:, i)' * N1 * x(:, i);
	term2(i) = K2 + x(:, i)' * L2 + x(:, i)' * N2 * x(:, i);
end
st_cov = cov(term1', term2'); % cov needs row vectors
st_cov = st_cov(1, 2); % only need the off diagonal element
an_cov = L1'*x_cov*L2 + 2*x_mean'*N2*x_cov*L1 + 2*x_mean'*N1*x_cov*L2 + 2*trace(N1*x_cov*N2*x_cov) + 4*x_mean'*N1*x_cov*N2*x_mean;
fprintf("Analytic: %f, Statistic: %f\n", an_cov, st_cov);