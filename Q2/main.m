%% Q1
clear; close all; clc;
mkdir results
addpath ../common/


% Constants
Kb = 1.38e-23;          % Boltzmann constant (J/K)
T0 = 290;               % Absolute temperature (K)
B = 1e6;                % Bandwidth (Hz)
F = 10^(5/10);          % Noise figure (linear scale)
L = 10^(3/10);          % System losses (linear scale)
SNR = 10^(10/10);       % Minimum SNR (linear scale)
lambda = 3e8 / 3.5e9;   % Wavelength (m), assuming frequency = 3.5 GHz
eta = 0.9;              % Efficiency factor
Nt = 4;                 % Number of transmitting antennas
Nr = 4;                 % Number of receiving antennas
Gt = 10^(10/10);        % Transmitting antenna gain (linear scale)
Gr = 10^(10/10);        % Receiving antenna gain (linear scale)

% Transmit power for mobile (W)
Pt_m = 0.0016;

% Calculate maximum LOS range R (m)
R = (lambda / (4 * pi)) * ...
    nthroot((eta * Nt * Pt_m * Nt * Gt * Nr * Gr) / ...
    (Kb * T0 * B * F * L * SNR), 2);

% Display result
fprintf('Maximum LOS Range: %.3f meters\n', R);

