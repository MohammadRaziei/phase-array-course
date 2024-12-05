function [fig1, fig2] = plot_af(theta, af)



fig1 = figure('units','normalized','outerposition',[0 .25 1 .5]); axis off
[subplot_axis, ~] = tight_subplot(1, 1, [0.2, 0.05], .1);
axes(subplot_axis(1)); axis on; hold on

plot(theta, af)
xlabel("\theta")
ylabel("Array Factor")
xlim([-1, 1]*pi)
xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))



af_db = mag2db(af); %20 log










fig2 = figure('units','normalized','outerposition',[0 .25 1 .5]); axis off
[subplot_axis, ~] = tight_subplot(1, 1, [0.2, 0.05], .1);
axes(subplot_axis(1)); axis on; hold on

plot(theta, af_db)
xlabel("\theta")
ylabel("Array Factor (db)")
xlim([-1, 1]*pi)
ylim([-60,60])
xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))

MaxIdx = find(islocalmax(af_db));

[af_db_max] = max(af_db);
hpbw_idx = find(af_db >= af_db_max - 3);
hpbw_idx = hpbw_idx((diff([0, hpbw_idx]) - diff([0 circshift(hpbw_idx, -1)]) ~= 0));

if hpbw_idx(1) == 1
    if length(hpbw_idx) > 4
        hpbw_idx = hpbw_idx(3:end-2);
    end
end


for i = MaxIdx
    text(theta(i), af_db(i)+15, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af_db(i),4))], "HorizontalAlignment", "left", "Color", .4*ones(1,3));
end

plot(theta(hpbw_idx), af_db(hpbw_idx), 'or')
plot(theta(MaxIdx), af_db(MaxIdx), '.r')

temp = af_db(MaxIdx);
peak = max(temp);
sll_value = max(temp(temp < .9 * peak));

sll = af_db_max - sll_value;
for i = hpbw_idx
    text(theta(i), af_db(i)+15, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af_db(i),4))], "HorizontalAlignment", "left", "Color", [1,.7,.7]);
end



for i = 1:length(hpbw_idx)/2
    drawbrace([theta(hpbw_idx(2*i)), af_db(hpbw_idx(2*i))-5], [theta(hpbw_idx(2*i-1)), af_db(hpbw_idx(1))-5], 0.01*diff(theta(hpbw_idx(2*i-1:2*i))), "Color", "red");
    HPBW = diff(theta(hpbw_idx(2*i-1:2*i)))/pi;
    text(mean(theta(hpbw_idx(2*i-1:2*i))), mean(af_db(hpbw_idx(2*i-1:2*i)))-20, sprintf("HPBW = %.2f\\pi", HPBW), "HorizontalAlignment", "center", "Color", [1,0,0]);
end


text(0, 45, [sprintf("Side Lobe level = %.2f db", sll); ...
    sprintf("peak = %.2g db", peak); ...
    sprintf("Half Power Band Width (HPBW) = %.2f\\pi", HPBW); ...
    ], "HorizontalAlignment", "center")
