clear; close all; clc;
mkdir results
addpath ../common/

w_n = [1, 1, 0, 0, 1, 0, 1];
d_lambda = 0.5;

AF_theta =@(w_n, d_lambda, theta_0) w_n * exp(1j*(0:length(w_n)-1).'* (2*pi*d_lambda * cos(theta_0)));
AF = @(w_n, d_lambda, theta) arrayfun(@(theta_0) AF_theta(w_n, d_lambda, theta_0), theta,'UniformOutput',true);

theta = -pi:0.001:pi;
af_a = abs(AF(w_n, d_lambda, theta));
af_b = abs(AF(ones(1,length(w_n)), d_lambda, theta));


figure; hold on; grid on
plot(theta, af_a, 'r')
plot(theta, af_b, 'k')
xlim([-1, 1]*pi);
xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))
exportgraphics(gcf, 'results/AF-plot.pdf', 'Append', false);
set(gca, 'YScale', 'log')
exportgraphics(gcf, 'results/AF-plot-logy.pdf', 'Append', false);


figure; 
polarplot(theta, af_b, 'k')
hold on;
polarplot(theta, af_a, 'r')
set(gca, 'ThetaZeroLocation', "top");
exportgraphics(gcf, 'results/AF-polarplot.pdf', 'Append', false);


%%
figure('units','normalized','outerposition',[0 .25 1 .5])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

grid on; hold on
xlim([-1, 1]*pi)
xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))
ylim([-2,8])
[~,MaxIdx] = findpeaks(af_a);
[~,MinIdx] = findpeaks(-af_a);

plot(theta, af_a, 'r')
plot(theta(MaxIdx), af_a(MaxIdx), 'xb')
plot(theta(MinIdx), af_a(MinIdx), 'ob')

for i = MaxIdx
    text(theta(i), af_a(i)+.5, sprintf("%.2f\\pi", theta(i)/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
end
for i = MinIdx
    text(theta(i), af_a(i)-.5, sprintf("%.2f\\pi", theta(i)/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
end

idx1 = MinIdx([floor(length(MinIdx) * 1 / 4), floor(length(MinIdx) * 1 / 4) + 1]);
idx2 = MinIdx([floor(length(MinIdx) * 3 / 4), floor(length(MinIdx) * 3 / 4) + 1]);
text(mean(theta(idx1)), mean(af_a(idx1))-1.4, sprintf("BW_{nn} = %.2f\\pi", diff(theta(idx1))/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3));
text(mean(theta(idx2)), mean(af_a(idx2))-1.4,  sprintf("BW_{nn} = %.2f\\pi", diff(theta(idx2))/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3));

drawbrace([theta(idx1(2)),af_a(idx1(2))-.8], [theta(idx1(1)),af_a(idx1(1))-.8], 0.006, "Color", "blue")
drawbrace([theta(idx2(2)),af_a(idx2(2))-.8], [theta(idx2(1)),af_a(idx2(1))-.8], 0.006, "Color", "blue")

exportgraphics(gcf, 'results/AF-plot-localmaxmin-a.pdf', 'Append', false);

%%
figure('units','normalized','position',[0 .25 1 .5])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];

grid on; hold on
xlim([-1, 1]*pi)
xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))
ylim([-2,8])
[~,MaxIdx] = findpeaks(af_b);
[~,MinIdx] = findpeaks(-af_b);

plot(theta, af_b, 'k')
plot(theta(MaxIdx), af_b(MaxIdx), 'xb')
plot(theta(MinIdx), af_b(MinIdx), 'ob')

for i = MaxIdx
    text(theta(i), af_b(i)+.5, sprintf("(%.2f\\pi, %g)", theta(i)/pi, round(af_b(i),2)), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
end
for i = MinIdx
    text(theta(i), af_b(i)-.5, sprintf("(%.2f\\pi, %g)", theta(i)/pi, round(af_b(i),2)), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
end
idx1 = MinIdx([floor(length(MinIdx) * 1 / 4), floor(length(MinIdx) * 1 / 4) + 1]);
idx2 = MinIdx([floor(length(MinIdx) * 3 / 4), floor(length(MinIdx) * 3 / 4) + 1]);
text(mean(theta(idx1)), mean(af_b(idx1))-1.4, sprintf("BW_{nn} = %.2f\\pi", diff(theta(idx1))/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3));
text(mean(theta(idx2)), mean(af_b(idx2))-1.4,  sprintf("BW_{nn} = %.2f\\pi", diff(theta(idx2))/pi), "HorizontalAlignment", "center", "Color", .3*ones(1,3));

drawbrace([theta(idx1(2)),af_b(idx1(2))-.8], [theta(idx1(1)),af_b(idx1(1))-.8], 0.006, "Color", "blue")
drawbrace([theta(idx2(2)),af_b(idx2(2))-.8], [theta(idx2(1)),af_b(idx2(1))-.8], 0.006, "Color", "blue")

exportgraphics(gcf, 'results/AF-plot-localmaxmin-b.pdf', 'Append', false);

%%



