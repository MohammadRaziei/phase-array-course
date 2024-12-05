%% Initialization
clear; close all; clc;
mkdir results
addpath ../common/ ../HW/



N = 21; % Number of points
d_lambda = 0.5;

AF_theta =@(w_n, d_lambda, theta_0) w_n * exp(1j*(0:length(w_n)-1).'* (2*pi*d_lambda * cos(theta_0)));
AF = @(w_n, d_lambda, theta) abs(arrayfun(@(theta_0) AF_theta(w_n, d_lambda, theta_0), theta,'UniformOutput',true));


% Generate Kaiser window
w = hamming_window(N);
csvwrite("results/hamming.csv", w');

% Display results
disp('Hamming Window:');
disp(w);

% Plot the window
figure;
stem(0:N-1, w, 'filled');
xlabel('Index (n)');
ylabel('Amplitude');
title('Hamming Window');
grid on;
ylim([0, 1.2])


exportgraphics(gcf, 'results/hamming-21.pdf', 'Append', false);


theta = -pi:0.001:pi;


af = AF(ones(1,N), d_lambda, theta);

[fig1, fig2] = plot_af(theta, af);

exportgraphics(fig1, 'results/af.pdf', 'Append', false);
exportgraphics(fig2, 'results/af-db.pdf', 'Append', false);


af = AF(w, d_lambda, theta);

[fig1, fig2] = plot_af(theta, af);
exportgraphics(fig1, 'results/af-with-hamming.pdf', 'Append', false);

exportgraphics(fig2, 'results/af-db-with-hamming.pdf', 'Append', false);

