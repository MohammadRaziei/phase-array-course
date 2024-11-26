function [fig1, fig2, fig3, theta_steering] = calc_beam_with_beta(beta, lambda0, lambda, is_plot)
if nargin == 1, lambda0_lambda = 1; else, lambda0_lambda = lambda0 / lambda; end
if nargin < 4, is_plot = true; end

N = 10;


AF_theta = @(N, lambda0_lambda, beta, theta_0) ones(1,N) * exp(1j*(0:N-1).'* (pi*lambda0_lambda * cos(theta_0)+beta));
AF = @(N, lambda0_lambda, beta, theta) arrayfun(@(theta_0) AF_theta(N, lambda0_lambda, beta, theta_0), theta,'UniformOutput',true);


theta = -pi:0.001:pi;


af = abs(AF(N, lambda0_lambda, beta, theta));

if is_plot
fig1 = figure;
polarplot(theta, af)
set(gca, "ThetaZeroLocation", "top");
title(sprintf("$\\Delta \\phi = %g\\,\\pi  \\ \\mathrm{rad}$", round(beta/pi,3)), "Interpreter", "latex");

exportgraphics(gcf, sprintf('results/AF-polarplot-beta-%i.pdf', round(rad2deg(beta))), 'Append', false);
else, fig1 = []; end


MaxIdx = find(islocalmax(af));

steering = [];
max_af = max(af);
i = 0;
for m = MaxIdx
    if af(m) > .9 * max_af
        i = i + 1;
        steering(i) = m;
    end
end

theta_steering = theta(steering);

if is_plot
fig2 = figure('units','normalized','outerposition',[0 .25 1 .5]); axis off
[subplot_axis, ~] = tight_subplot(1, 1, [0.2, 0.05], .1);
axes(subplot_axis(1)); axis on; hold on

plot(theta, af)


plot(theta(MaxIdx), af(MaxIdx), 'xb')

for i = MaxIdx
    text(theta(i), af(i)+1, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af(i),4))], "HorizontalAlignment", "left", "Color", [.4,.4,.4]);
end

for t = theta_steering
    xline(t);
    text(t, 1, [sprintf("$\\theta = %g\\pi$", t/pi); sprintf("$\\theta = %g^\\circ$", round(rad2deg(t),1))], "Interpreter", "latex", "HorizontalAlignment", "left");
end

x_range_pi = -1:0.125:1;
xlim([-pi, pi]);
xticks(x_range_pi*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), x_range_pi,'UniformOutput',true))
xlabel("$\theta$", "Interpreter", "latex")
title(sprintf("$\\Delta \\phi = %g\\,\\pi  \\ \\mathrm{rad}$", round(beta/pi,3)), "Interpreter", "latex");
ylim([0, 12])


exportgraphics(gcf, sprintf('results/AF-plot-beta-%i.pdf', round(rad2deg(beta))), 'Append', false);
else, fig2 = []; end
%%
af_db = mag2db(af);


MaxIdx = find(islocalmax(af_db));


[af_db_max] = max(af_db);
hpbw_idx = find(af_db >= af_db_max - 3);
hpbw_idx = hpbw_idx((diff([0, hpbw_idx]) - diff([0 circshift(hpbw_idx, -1)]) ~= 0));


if hpbw_idx(1) == 1
    if length(hpbw_idx) > 4
        hpbw_idx = hpbw_idx(3:end-2);
    end
end

if is_plot
fig3 = figure('units','normalized','outerposition',[0 .25 1 .5]); axis off
[subplot_axis, ~] = tight_subplot(1, 1, [0.2, 0.05], .1);
axes(subplot_axis(1)); axis on; hold on

plot(theta, af_db)


plot(theta(MaxIdx), af_db(MaxIdx), 'xb')

for i = MaxIdx
    text(theta(i), af_db(i)+15, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af_db(i),4))], "HorizontalAlignment", "left", "Color", .4*ones(1,3));
end


plot(theta(hpbw_idx), af_db(hpbw_idx), 'or')


for i = hpbw_idx
    text(theta(i), af_db(i)+15, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af_db(i),4))], "HorizontalAlignment", "left", "Color", [1,.7,.7]);
end

for i = 1:length(hpbw_idx)/2
    drawbrace([theta(hpbw_idx(2*i)), af_db(hpbw_idx(2*i))-5], [theta(hpbw_idx(2*i-1)), af_db(hpbw_idx(1))-5], 0.01*diff(theta(hpbw_idx(2*i-1:2*i))), "Color", "red");
    text(mean(theta(hpbw_idx(2*i-1:2*i))), mean(af_db(hpbw_idx(2*i-1:2*i)))-20, sprintf("HPBW = %.2f\\pi", diff(theta(hpbw_idx(2*i-1:2*i)))/pi), "HorizontalAlignment", "center", "Color", [1,0,0]);
end

xlim([-pi, pi]);
xticks(x_range_pi*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), x_range_pi,'UniformOutput',true))
xlabel("$\theta$", "Interpreter", "latex")
ylim([-100, 70])
title(sprintf("$\\Delta \\phi = %g\\,\\pi  \\ \\mathrm{rad}$", round(beta/pi,3)), "Interpreter", "latex");


exportgraphics(gcf, sprintf('results/AF-plot-db-beta-%i.pdf', round(rad2deg(beta))), 'Append', false);
else, fig3 = []; end

