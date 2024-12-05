%% Q4
clear; close all; clc;
mkdir results
addpath ../common/
addpath ../HW/

rng(1)


AF = @(P, lambda, theta, phi) abs(cell2mat(arrayfun(@(phi) arrayfun(@(theta) sum(exp(1j*2*pi/lambda *[sin(theta)* cos(phi), sin(theta)* sin(phi), cos(theta)]*P')), theta,'UniformOutput',true)', phi,'UniformOutput',false)));


theta_d = deg2rad(30);
phi_d = deg2rad(45);

%%

lambda = 1;
Nx = 4;
Ny = 4;
dx = .6;
dy = .6;
P = cell2mat(arrayfun(@(x) [(floor(x/Ny)-(Nx-1)/2)*dx, (mod(x,Ny)-(Ny-1)/2)*dy, 0]', 0:Nx*Ny-1,'UniformOutput',false))';

weights = ones(Nx, Ny);

phi_s = phi_d;

plot_beamforming(P, weights, phi_s, "c");

%%


% Initialize the matrix
weights = zeros(Nx, Ny);
% Loop to calculate the weights
for m = 0:Nx-1
    for n = 0:Ny-1
        weights(m+1, n+1) = exp(-2j * pi * (m * sin(theta_d)*cos(phi_d) * dx / lambda + ...
            n * sin(theta_d)*sin(phi_d) * dy / lambda));
    end
end

% Display the matrix
disp('Matrix of weights (W):');
disp(weights);

phi_s = phi_d;

plot_beamforming(P, weights, phi_s, "w");

%%

lambda = 1;
Nx = 5;
Ny = 3;
dx = .5;
dy = .8;
P = cell2mat(arrayfun(@(x) [(floor(x/Ny)-(Nx-1)/2)*dx, (mod(x,Ny)-(Ny-1)/2)*dy, 0]', 0:Nx*Ny-1,'UniformOutput',false))';

% Initialize the matrix
weights = zeros(Nx, Ny);
% Loop to calculate the weights
for m = 0:Nx-1
    for n = 0:Ny-1
        weights(m+1, n+1) = exp(-2j * pi * (m * sin(theta_d)*cos(phi_d) * dx / lambda + ...
            n * sin(theta_d)*sin(phi_d) * dy / lambda));
    end
end

% Display the matrix
disp('Matrix of weights (W):');
disp(weights);

phi_s = phi_d;

plot_beamforming(P, weights, phi_s, "w2");

