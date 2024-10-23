clear; close all; clc;
mkdir results
addpath ../common/
addpath ../HW/


R1=(2.5e-3)*...
[3 0 0
 2 0 0
 1 0 0
 0 0 0
-1 0 0
-2 0 0
-3 0 0];



phi = 0;
phases = zeros(1, size(R1, 1));
fs = 100;
theta = -180:1/fs:180;
arrayfactor = arrayfun(@(theta) Array_beam_cal(R1,phases,theta,phi), theta,'UniformOutput',true);
theta = deg2rad(theta);

figure
polarplot(theta, arrayfactor)

figure; hold on
plot(deg2rad(theta), arrayfactor)


af_smooth = smooth(arrayfactor, 4*fs);
plot(deg2rad(theta), af_smooth)


[~,MaxIdx] = findpeaks(af_smooth);
[~,MinIdx] = findpeaks(-af_smooth);

% plot(theta, arrayfactor, 'r')
% plot(theta(MaxIdx), arrayfactor(MaxIdx), 'xb')
% plot(theta(MinIdx), arrayfactor(MinIdx), 'ob')
% 
% for i = MaxIdx
%     text(theta(i), arrayfactor(i)+.5, sprintf("%.2f\\pi", theta(i)/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
% end
% for i = MinIdx
%     text(theta(i), arrayfactor(i)-.5, sprintf("%.2f\\pi", theta(i)/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
% end

MaxIdx



R3=(5e-3)*...
[3+0.5*randn() 0 0
2+0.5*randn() 0 0
1+0.5*randn() 0 0
0+0.5*randn() 0 0
-1+0.5*+randn() 0 0
-2+0.5*randn() 0 0
-3+0.5*randn() 0 0];



R2=(5e-3)*...
[3 0 0
2 0 0
1 0 0
30 0 0
-1 0 0
-2 0 0
-3 0 0];