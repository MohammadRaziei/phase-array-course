%% Q1
clear; close all; clc;
mkdir results
addpath ../common/
% Initialize parameters


M = 4;                 % Modulation order
bps = log2(M);         % Bits per symbol
N = 7;                 % RS codeword length
K = 5;                 % RS message length

% Define simulation parameters
ebnoVec = (0:0.006:15)'; % Range of Eb/N0 values in dB
targetBER = 1e-5;      % Target BER threshold

% Calculate BER with Reed-Solomon coding
berwCoding = bercoding(ebnoVec, 'RS', 'hard', N, K, 'psk', M, 'nondiff');

% Find the index of the last BER value above the threshold
idx = max(find(berwCoding >= targetBER));

% Plot results
figure;
hold on;
set(gca, 'YScale', 'log');
plot(ebnoVec, berwCoding, 'b', 'LineWidth', 1.5); % BER curve
plot(ebnoVec(idx), berwCoding(idx), 'ro', 'MarkerSize', 5); % Threshold point
yline(targetBER, 'r--', 'LineWidth', 1); % Target BER line

% Annotate the threshold point
text(ebnoVec(idx), berwCoding(idx)*10, ...
     sprintf('(%g, %g)', ebnoVec(idx), berwCoding(idx)), ...
     'FontSize', 10, 'HorizontalAlignment', 'left');

% Labels and title
xlabel('SNR or Eb/N0 (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR with RS Coding');
grid on;

% Export the plot
exportgraphics(gcf, 'results/ber-snr.pdf', 'Append', false);

%%
saver = [];
snr_db = ebnoVec(idx);
saver(1,:) = snr_db;

snr = db2pow(snr_db);
saver(2,:) = snr;

R = 1e9;
C = 1.25*R;
Bw = C / log2(1+snr);
saver(3,:) = Bw;

d = 1.2e3;
L_db = 18;
L = db2mag(L_db);
saver(4,:) = L;

f = 60e9;
lambda = 3e8 / f;
saver(5,:) = lambda;

Kb = 1.38e-23;
T0 = 290;

F = 5;
Pt = 8;

G2 = Kb * T0 * Bw * snr * L * db2pow(F) / db2pow(Pt) * (4 * pi * d / lambda)^2;
saver(6,:) = G2;
G = sqrt(G2);
saver(7,:) = G;
G_db = pow2db(G);
saver(8,:) = G_db;






csvwrite("results/saver.csv", saver);











