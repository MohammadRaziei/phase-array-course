function run_for_sll(sll)

w = dolph_chebyshev(11, sll);

csvwrite(sprintf("results/w-sll-%g.csv", sll), round(w',4));

AF = @(w_n, theta) arrayfun(@(theta_0) w_n * exp(1j*(0:length(w_n)-1).'* (pi* cos(theta_0))), theta,'UniformOutput',true);

theta = -pi:0.001:pi;

af = abs(AF(w, theta));

MaxIdx = find(islocalmax(af));

figure
polarplot(theta, af);
exportgraphics(gcf, sprintf("results/af-polar-sll-%g.pdf", sll), 'Append', false);

af_db = mag2db(af);

figure('units','normalized','outerposition',[0 .25 1 .5]); axis off
[subplot_axis, ~] = tight_subplot(1, 1, [0.2, 0.05], .1);
axes(subplot_axis(1)); axis on; hold on

plot(theta, af_db);

x_range_pi = -1:0.125:1;
xlim([-pi, pi]);
xticks(x_range_pi*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), x_range_pi,'UniformOutput', true))
xlabel("$\theta$", "Interpreter", "latex")
ylim([-60, 50])

for i = MaxIdx
    text(theta(i), af_db(i)+10, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af_db(i),4))], "HorizontalAlignment", "left", "Color", [.4,.4,.4]);
end

[af_db_max] = max(af_db);
hpbw_idx = find(af_db >= af_db_max - 3);
hpbw_idx = hpbw_idx((diff([0, hpbw_idx]) - diff([0 circshift(hpbw_idx, -1)]) ~= 0));


if hpbw_idx(1) == 1
    if length(hpbw_idx) > 4
        hpbw_idx = hpbw_idx(3:end-2);
    end
end

plot(theta(hpbw_idx), af_db(hpbw_idx), 'or')

for i = hpbw_idx
    text(theta(i), af_db(i)+8, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af_db(i),4))], "HorizontalAlignment", "left", "Color", [1,.7,.7]);
end

for i = 1:length(hpbw_idx)/2
    drawbrace([theta(hpbw_idx(2*i)), af_db(hpbw_idx(2*i))-3], [theta(hpbw_idx(2*i-1)), af_db(hpbw_idx(1))-3], 0.01*diff(theta(hpbw_idx(2*i-1:2*i))), "Color", "red");
    text(mean(theta(hpbw_idx(2*i-1:2*i))), mean(af_db(hpbw_idx(2*i-1:2*i)))-10, sprintf("HPBW = %.2f\\pi", diff(theta(hpbw_idx(2*i-1:2*i)))/pi), "HorizontalAlignment", "center", "Color", [1,0,0]);
end

MinIdx = find(islocalmin(af));
for i = MinIdx
    text(theta(i), mean([-20, af_db(i)]), sprintf("%.2f\\pi", theta(i)/pi), "HorizontalAlignment", "center", "Color", [.4,.4,.4]);
end

exportgraphics(gcf, sprintf("results/af-sll-%g.pdf", sll), 'Append', false);



