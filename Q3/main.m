%% Initialization
clear; close all; clc;
mkdir results
addpath ../common/ ../HW/

%% Load Data
data = readtable('Tpatch.csv');
phi_angles = data.Phi_deg_;
theta_angles = data.Theta_deg_;
gain_total_db = data.dB_GainTotal_;
pattern_reshaped = reshape(gain_total_db, 181, []); % Reshape for 181 x 181 grid
theta_grid = (min(theta_angles):1:max(theta_angles)).';
phi_grid = (min(phi_angles):1:max(phi_angles)).';

%% Constants
c = 3e8;              % Speed of light (m/s)
Kb = 1.38e-23;        % Boltzmann constant (J/K)
T0 = 290;             % Absolute temperature (K)
f0 = 28e9;            % Operating frequency (Hz)
lambda = c / f0;      % Wavelength (m)
element_spacing = lambda / 2; % Element spacing (m)
Gt_db = 18;           % Transmit antenna gain (dB)
Gt = 10^(Gt_db / 10); % Transmit antenna gain (linear scale)
Gr_db = 18;           % Receive antenna gain (dB)
Gr = 10^(Gr_db / 10); % Receive antenna gain (linear scale)

%% Array Factor Calculations
N_elements = 4; % Number of elements in x-direction
M_elements = 8; % Number of elements in z-direction

% Compute Array Factor in x-direction
sincos_component = cosd(phi_grid) * sind(theta_grid.');
AF_x = sin(N_elements / 2 * 2 * pi / lambda * element_spacing * sincos_component) ./ ...
       sin(0.5 * 2 * pi / lambda * element_spacing * sincos_component);
AF_x(isnan(AF_x)) = N_elements; % Handle NaN cases

% Compute Array Factor in z-direction
theta_repeated = repmat(theta_grid, 1, length(phi_grid));
AF_z = sin(M_elements / 2 * 2 * pi / lambda * element_spacing * cosd(theta_repeated')) ./ ...
       sin(0.5 * 2 * pi / lambda * element_spacing * cosd(theta_repeated'));
AF_z(isnan(AF_z)) = M_elements; % Handle NaN cases

% Total Array Factor
array_factor_db = 10 * log10(abs(AF_x .* AF_z));

%% Total Pattern
total_pattern_db = array_factor_db + pattern_reshaped;
plotPattern(10.^(total_pattern_db / 10), theta_grid, phi_grid); % Plot pattern
title("\textbf{Pattern of 4 $\times$ 8 Planner Array}", 'Interpreter', 'latex');
G_max = max(total_pattern_db, [], 'all'); % Maximum gain

% Export to PDF
exportgraphics(gcf, 'results/Pattern-Plot-4x8.pdf', 'Append', false);


%% Link Budget Calculations
eta = 1;             % Efficiency factor
L_db = 0.04;         % System loss (dB)
L = 10^(L_db / 10);  % System loss (linear scale)
F_db = 10;           % Noise figure (dB)
F = 10^(F_db / 10);  % Noise figure (linear scale)
SNR = 10;            % Signal-to-noise ratio (linear scale)
B = 361e6;           % Bandwidth (Hz)
distance = 200;      % Link distance (m)
n = 2;               % Path loss exponent
Gt_db = 22.42;       % Updated transmit antenna gain (dB)
Gt = 10^(Gt_db / 10);% Updated transmit antenna gain (linear scale)
Gr_db = 19.41;       % Updated receive antenna gain (dB)
Gr = 10^(Gr_db / 10);% Updated receive antenna gain (linear scale)
Nt = 32;             % Number of transmitting antennas
Nr = 16;             % Number of receiving antennas

% Transmit Power
Pt = SNR * Kb * T0 * B * F * L / (eta * Gt * Gr * Nt) * (4 * pi * distance / lambda)^n;
Pt_dbm = 10 * log10(Pt * 1000); % Convert to dBm
fprintf('Transmit Power: %.3f W (%.3f dBm)\n', Pt, Pt_dbm);