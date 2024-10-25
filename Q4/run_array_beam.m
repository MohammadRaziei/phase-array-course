function [arrayfactor, fig1, fig2, MinIdx, MaxIdx] = run_array_beam(R1, N, is_plot)
if nargin < 3, is_plot = true; end
phi = 0;
phases = zeros(1, size(R1, 1));
fs = 100;
theta = -180:1/fs:180;
arrayfactor = arrayfun(@(theta) Array_beam_cal(R1,phases,theta,phi), theta,'UniformOutput',true);

theta = deg2rad(theta);

if is_plot, fig1 = figure;
polarplot(theta, arrayfactor)
set(gca, 'ThetaZeroLocation', "top");
else,fig1 = []; 
end

af_smooth = [arrayfactor, arrayfactor, arrayfactor];
for i = 1:N, af_smooth = smooth(af_smooth, 4*fs); end
af_smooth = af_smooth(length(arrayfactor):length(arrayfactor)*2-1);





[~,MaxIdx] = findpeaks(af_smooth);
[~,MinIdx] = findpeaks(-af_smooth);
MinIdx_ = MinIdx; MaxIdx_ = MaxIdx;
win = 5*fs;
for i = 1:length(MinIdx)
    midx = MinIdx(i);
    midx = (1:win)-floor(win/2)+midx;
    midx = midx(midx>0);
    [~,idx] = min(arrayfactor(midx));
    MinIdx(i) = midx(idx);
end
MinIdx = unique(MinIdx).';
for i = 1:length(MaxIdx)
    midx = MaxIdx(i);
    midx = (1:win)-floor(win/2)+midx;
    midx = midx(midx>0);
    [~,idx] = max(arrayfactor(midx));
    MaxIdx(i) = midx(idx);
end
MaxIdx = unique(MaxIdx).';

if is_plot, fig2 = figure('units','normalized','outerposition',[0 .25 1 .5]);
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + 2.2*ti(1);
bottom = outerpos(2) + 2.2*ti(2);
ax_width = outerpos(3) - 2.2*ti(1) - ti(3);
ax_height = outerpos(4) - 2.2*ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
hold on

plot(theta, arrayfactor)
plot(theta, af_smooth, "Color", ones(1, 3)*.5);
plot(theta(MaxIdx), arrayfactor(MaxIdx), 'xr')
plot(theta(MinIdx), arrayfactor(MinIdx), 'or')
plot(theta(MinIdx_), af_smooth(MinIdx_), '.b')
plot(theta(MaxIdx_), af_smooth(MaxIdx_), '.b')
ylim([-2, 10])
xlim([-pi, pi])
xticks((-1:0.25:1)*pi)
xlabel("\theta (rad)"); ylabel("Pattern (magnitude)")
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))

for i = MaxIdx
    text(theta(i), arrayfactor(i)+.9, sprintf("\\lceil%.2g\\pi", round(theta(i)/pi, 2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
    text(theta(i), arrayfactor(i)+.5, sprintf("\\lfloor%.2g", round(arrayfactor(i),2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
end
for i = MinIdx
    text(theta(i), arrayfactor(i)-.5, sprintf("\\lceil%.2g\\pi", round(theta(i)/pi, 2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
    text(theta(i), arrayfactor(i)-.9, sprintf("\\lfloor%.2g", round(arrayfactor(i),2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
end

else, fig2 = [];
end

