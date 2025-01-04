%% Q1
clear; close all; clc;
mkdir results
addpath ../common/ ../HW/
% Initialize parameters

run VGAS21

voltage = S21(:,1);
amplitude = S21(:,2);
phases = S21(:,3);

figure
plot(voltage, amplitude, "LineWidth", 1, "Color", "red")
xlabel("Control Voltage"); ylabel("Gain Amplitude")
ylim([0,4])
grid on
exportgraphics(gca, "results/amplitude-vs-voltage.pdf", "Append", false);

figure
plot(voltage, phases, "LineWidth", 1, "Color", "blue")
xlabel("Control Voltage"); ylabel("Phase (deg)")
grid on
exportgraphics(gca, "results/phase-vs-voltage.pdf", "Append", false);



% Vc = VcI + 1i * VcQ;
%%
% dG = min(abs(diff(amplitude)));
indices = round(linspace(1, numel(amplitude), 8)); % Get indices of 8 elements
% idx = numel(amplitude)-8+1:numel(amplitude);
g = amplitude(indices);

gI = g;
gQ = g;
[gI, gQ] = meshgrid(gI, gQ); 
gain = gI + 1i * gQ;

figure
hold on
xlim([0,3.2]); ylim([0, 3.2])

scatter(real(gain), imag(gain), '.', "MarkerEdgeColor", "black");
xlabel("Vc_I"); ylabel("Vc_Q")


GainFlat = gain(:);
a = abs(GainFlat);
Rg = 3;
dG = .52;

gI = 0:.02:1.5*max(gI(:));
plot(gI, sqrt((Rg+dG)^2 - (gI).^2), ":b");
plot(gI, sqrt((Rg-dG)^2 - (gI).^2), ":b");


idx = find(and(a < Rg + dG, a > Rg - dG));

GainIdx = GainFlat(idx);
n = numel(GainIdx);

plot(real(GainIdx), imag(GainIdx), "or")
plot([real(GainIdx) zeros(n,1)]', [imag(GainIdx) zeros(n,1)]', ":k")

exportgraphics(gca, "results/vga-phase-map.pdf", "Append", false);

disp(n)

figure; hold on
phase = (sort(angle(GainIdx)));
phase = rad2deg(phase);
plot(phase);
plot(phase, ".");
ylabel("Phase (deg)")
exportgraphics(gca, "results/vga-phase-sort.pdf", "Append", false);





%%
indices = round(linspace(1, numel(amplitude), 8)); % Get indices of 8 elements
g = amplitude(indices) .* exp(1i * deg2rad(phases(indices)));

gI = g;
gQ = g;
[gI, gQ] = meshgrid(gI, gQ); 
gain = gI + 1i * gQ;

figure
hold on
xlim([-3.2, 1.5]); ylim([0, 3.7])

scatter(real(gain), imag(gain), '.', "MarkerEdgeColor", "black");
xlabel("Vc_I"); ylabel("Vc_Q")


GainFlat = gain(:);
a = abs(GainFlat);
Rg = 3;
dG = .52;

gI = -3.2:.02:3;
plot(gI, sqrt((Rg+dG)^2 - (gI).^2), ":b");
plot(gI, sqrt((Rg-dG)^2 - (gI).^2), ":b");


idx = find(and(a < Rg + dG, a > Rg - dG));

GainIdx = GainFlat(idx);
n = numel(GainIdx);

plot(real(GainIdx), imag(GainIdx), "or")
plot([real(GainIdx) zeros(n,1)]', [imag(GainIdx) zeros(n,1)]', ":k")

exportgraphics(gca, "results/vga-phase-map-nonideal.pdf", "Append", false);

disp(n)

figure; hold on
phase = (sort(angle(GainIdx)));
phase = rad2deg(phase);
plot(phase);
plot(phase, ".");
ylabel("Phase (deg)")
exportgraphics(gca, "results/vga-phase-sort-nonideal.pdf", "Append", false);





