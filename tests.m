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
	term = term + (x(:, i)' * B * x(:, i)) * (x(:, i)' * A);
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