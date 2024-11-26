%% Q5
clear; close all; clc;
mkdir results
addpath ../common/


d_lambda = .5;
N = 10;


AF_theta =@(N, d_lambda, beta, theta_0) ones(1,N) * exp(1j*(0:N-1).'* (2*pi*d_lambda * cos(theta_0)+beta));
AF = @(N, d_lambda, beta, theta) arrayfun(@(theta_0) AF_theta(N, d_lambda, beta, theta_0), theta,'UniformOutput',true);


theta = -pi:0.001:pi;


%%
beta = 0;
calc_beam_with_beta(beta);

%%
beta = pi/2;
calc_beam_with_beta(beta);

%%
beta = pi;
calc_beam_with_beta(beta);


%%
beta = pi/4;
c = 3e8;

f0 = 60e9;
lambda0 = c / f0;
f = 60e9;
lambda = c / f;
calc_beam_with_beta(beta, lambda0, lambda, true);
%%
freqs = linspace(50e9, 70e9, 100);
theta_steerings = zeros(numel(freqs), 2);
for i = 1:numel(freqs)
    lambda = c / freqs(i);
    [~, ~, ~, t] = calc_beam_with_beta(beta, lambda0, lambda, false);
    theta_steerings(i,:) = rad2deg(t);
end

figure
yyaxis left
plot(freqs/1e9, theta_steerings(:,1))
ylabel("$\theta$ (degree)", "Interpreter", "latex")

yyaxis right
plot(freqs/1e9, theta_steerings(:,2))
ylabel("$\theta$ (degree)", "Interpreter", "latex")

xlabel("frequency (GHz)", "Interpreter", "latex")

title("Streeing angle for $\frac{\pi}{4}\ \mathrm{rad}$ phase shift for $f_0 = 60\ \mathrm{GHz}$", "Interpreter","latex");

exportgraphics(gcf, "results/streeing-45.pdf", 'Append', false);




%%
%%
freqs = linspace(50e9, 70e9, 100);
theta_steerings = zeros(numel(freqs), 2);
for i = 1:numel(freqs)
    lambda = c / freqs(i);
    [~, ~, ~, t] = calc_beam_with_beta(2*pi*lambda0/lambda, lambda0, lambda, false);
    theta_steerings(i,:) = rad2deg(t);
end

figure
yyaxis left
plot(freqs/1e9, theta_steerings(:,1))
ylabel("$\theta$ (degree)", "Interpreter", "latex")

yyaxis right
plot(freqs/1e9, theta_steerings(:,2))
ylabel("$\theta$ (degree)", "Interpreter", "latex")

xlabel("frequency (GHz)", "Interpreter", "latex")

title("Streeing angle for $\frac{\pi}{4}\ \mathrm{rad}$ phase shift for $f_0 = 60\ \mathrm{GHz}$", "Interpreter","latex");

exportgraphics(gcf, "results/streeing-d.pdf", 'Append', false);


