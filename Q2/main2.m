%% Q4
clear; close all; clc;
mkdir results
addpath ../common/
addpath ../HW/

rng(1)


AF = @(P, lambda, theta, phi) abs(cell2mat(arrayfun(@(phi) arrayfun(@(theta) sum(exp(1j*2*pi/lambda *[sin(theta)* cos(phi), sin(theta)* sin(phi), cos(theta)]*P')), theta,'UniformOutput',true)', phi,'UniformOutput',false)));



N = 21;
lambda = 1;
d = lambda/2;


%% a

P = [zeros(N, 1), zeros(N, 1), (0:N-1)'*d];


run_for_P(P, "a");



%% b

x = linspace(0, 2*pi, N+1)';
[x, y, z] = sph2cart(x(1:end-1), zeros(N, 1), ones(N, 1));
P = [x, y, z];
% Compute the distance matrix
D = pdist(P);           % Pairwise distances as a vector
min_dist = min(D);      % Find the minimum distance
P = P * d / min_dist;
clear x y z D min_dist

run_for_P(P, "b");






