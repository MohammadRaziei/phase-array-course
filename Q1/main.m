%% Q1
clear; close all; clc;
mkdir results
addpath ../common/ 
% Initialize parameters

AmpW = load("../HW/AmpW.mat");
PhasW = load("../HW/PhasW.mat");

% Freq, Vcontrol, Phaseshift

PhaseShifterIdx = [2, 3, 8, 17];
vControl = 0:10;

%%
% Determine the middle index along the frequency axis
midFreqIdx = ceil(size(AmpW.PV2, 1) / 2);

% Extract amplitude values for the specific frequency index
amplitudes = squeeze(AmpW.PV2(midFreqIdx, :, PhaseShifterIdx));

% Normalize amplitudes by the maximum value for each phase shift
normalizedAmplitudes = amplitudes ./ max(amplitudes, [], 1);

% Plot the normalized amplitudes vs Vcontrol



figure('units','normalized','position',[0 .25 1 .35]); hold on
plot(normalizedAmplitudes, '-o', 'LineWidth', 1.5, 'MarkerSize',3);
xlabel('Vcontrol');
ylabel('Normalized Amplitude');
title('Phase Shifter Amplitude vs Vcontrol');
legend(arrayfun(@(i) sprintf('Phase Shift %d', i), 1:size(amplitudes, 2), 'UniformOutput', false), 'Location', 'bestoutside'); % Update based on number of phase shifts
grid on;
exportgraphics(gca, "results/amplitudeVcontrol.pdf", "Append", false)


% Determine the middle index along the frequency axis
midFreqIdx = ceil(size(PhasW.PV2, 1) / 2);

% Extract phase values for the specific frequency index
phases = squeeze(PhasW.PV2(midFreqIdx, :, PhaseShifterIdx));
% Plot the phase values vs Vcontrol
figure('units','normalized','position',[0 .25 1 .35]); hold on;
plot(vControl, phases, '-o', 'LineWidth', 1.5, 'MarkerSize',3, 'MarkerFaceColor','auto');
xlabel('Vcontrol');
ylabel('Phase Shift (radians)');
title('Phase Shift vs Vcontrol');
legend(arrayfun(@(i) sprintf('Phase Shift %d', i), 1:size(amplitudes, 2), 'UniformOutput', false), 'Location', 'bestoutside'); % Update based on number of phase shifts
grid on;
exportgraphics(gca, "results/phaseVcontrol.pdf", "Append", false)








%%
% Define PhaseShifterIdx and control arrays
PhaseShifterIdx = [2, 3, 8, 17];
vControlOrig = 0:10; 
vControl = linspace(0, 10, 100);
vControl4bit = linspace(0, 10, 2^4); % 4-bit control levels

%% Amplitude Plot
% Interpolation using interp1 with 'cubic' for amplitudes
interpAmplitudes = interp1(vControlOrig, normalizedAmplitudes, vControl, 'cubic');
amplitudes4bit = interp1(vControlOrig, normalizedAmplitudes, vControl4bit, 'cubic');

% Plot normalized amplitudes and interpolated values
figure('units','normalized','position',[0 .25 1 .35]); hold on;
plot(vControl, interpAmplitudes, '-', 'LineWidth', 1.5); % Original points
plot(vControl4bit, amplitudes4bit, 'k.', 'MarkerSize',8); % Interpolated points
xlabel('Vcontrol');
ylabel('Normalized Amplitude');
title('Phase Shifter Amplitude vs Vcontrol');
legend(arrayfun(@(i) sprintf('Phase Shift %d', i), 1:size(amplitudes, 2), 'UniformOutput', false), 'Location', 'bestoutside'); % Update based on number of phase shifts
grid on;
exportgraphics(gca, "results/amplitudeVcontrol4bit.pdf", "Append", false)



%% Phase Plot
% Interpolation using interp1 with 'cubic' for phases
interpPhases = interp1(vControlOrig, phases, vControl, 'cubic');
phases4bit = interp1(vControlOrig, phases, vControl4bit, 'cubic');

% Plot phase values and interpolated values
figure('units','normalized','position',[0 .25 1 .35]); hold on;
plot(vControl, interpPhases, '-', 'LineWidth', 1.5); % Interpolated curve
plot(vControl4bit, phases4bit, 'k.', 'MarkerSize', 8); % Highlight 4-bit control levels
xlabel('Vcontrol');
ylabel('Phase Shift (radians)');
title('Phase Shift vs Vcontrol');
legend(arrayfun(@(i) sprintf('Phase Shift %d', i), 1:size(phases, 2), 'UniformOutput', false), 'Location', 'bestoutside'); % Update based on number of phase shifts
grid on;
exportgraphics(gca, "results/phaseVcontrol4bit.pdf", "Append", false)


%%






c = 3e8;        % signal propagation speed
fc = 60e9;       % signal carrier frequency
lambda = c/fc;  % wavelength

thetaad = 3;     % look directions
thetaan = -90:5:90;           % interference direction

ula = phased.ULA(numel(PhaseShifterIdx), lambda/2);
ula.Element.BackBaffled = true;

D = getElementPosition(ula)/lambda;

weights = cell2mat(arrayfun(@(theta) nullstreeing(D, theta, thetaad).', thetaan.', "UniformOutput", false));



% Plot the pattern
% pattern(ula,fc,-180:180,0,'PropagationSpeed',c,'Type','powerdb',...
%     'CoordinateSystem','rectangular','Weights',weights.');
% hold on; legend off;
% plot([40 40],[-100 0],'r--','LineWidth',2)
% text(40.5,-5,'\leftarrow Interference Direction','Interpreter','tex',...
%     'Color','r','FontSize',10)

mkdir results/coherent

figure('units','normalized','position',[0 .25 1 .5]); 
for i = 1:size(weights)
    w = weights(i,:).';
    pattern(ula,fc,-180:180,0,'PropagationSpeed',c,'Type','powerdb',...
        'CoordinateSystem','rectangular','Weights',w);
    hold on; legend off;
    plot([thetaad thetaad],[-100 10],'g--','LineWidth',2)
    plot([thetaan(i) thetaan(i)],[-100 10],'r--','LineWidth',2)
    if thetaan(i) <= 60
        text(thetaan(i)+.5,5,'\leftarrow Interference Direction','Interpreter','tex',...
            'Color','r','FontSize',10)
    else
        text(thetaan(i)-.5,5,'Interference Direction \rightarrow','Interpreter','tex',...
            'HorizontalAlignment', 'right', 'Color','r','FontSize',10)
    end
    xlim([-91,91])
    exportgraphics(gca, sprintf("results/coherent/frame-%i.pdf", i), "Append", false)
    clf
end

fileID = fopen('results/coherent/smallsize.tml', 'w');
fprintf(fileID, '::0x0,1\n');
for d = 2:i
    fprintf(fileID, '::%d\n', d);
end
fclose(fileID);




% % Zoom
% xlim([30 50])
% legend(arrayfun(@(k)sprintf('%d degrees',k),thetaad,...
%     'UniformOutput',false),'Location','SouthEast')

