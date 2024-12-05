%% Initialization
clear; close all; clc;
mkdir results
addpath ../common/ ../HW/



N = 11; % Number of points
alpha = 3.427; % Alpha value
d_lambda = 0.5;

AF_theta =@(w_n, d_lambda, theta_0) w_n * exp(1j*(0:length(w_n)-1).'* (2*pi*d_lambda * cos(theta_0)));
AF = @(w_n, d_lambda, theta) abs(arrayfun(@(theta_0) AF_theta(w_n, d_lambda, theta_0), theta,'UniformOutput',true));


% Generate Kaiser window
w = kaiser_window(N, alpha);

% Display results
disp('Kaiser Window:');
disp(w);

% Plot the window
figure;
stem(0:N-1, w, 'filled');
xlabel('Index (n)');
ylabel('Amplitude');
title('Kaiser Window');
grid on;
ylim([0, 1.2])


exportgraphics(gcf, 'results/keiser-11.pdf', 'Append', false);


theta = -pi:0.001:pi;

af = AF(w, d_lambda, theta);


figure('units','normalized','outerposition',[0 .25 1 .5]); axis off
[subplot_axis, ~] = tight_subplot(1, 1, [0.2, 0.05], .1);
axes(subplot_axis(1)); axis on; hold on

plot(theta, af)
xlabel("\theta")
ylabel("Array Factor")
xlim([-1, 1]*pi)
xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))

exportgraphics(gcf, 'results/af-with-keiser.pdf', 'Append', false);
csvwrite("results/keiser.csv", w');


af_db = mag2db(af); %20 log


figure('units','normalized','outerposition',[0 .25 1 .5]); axis off
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

sll_value = -inf;
for i = 1:MaxIdx
    if af_db(i) < .9* af_db_max
        temp = af_db(i);
        if sll_value < temp
            sll_value = temp;
        end
    end
end
temp = af_db(MaxIdx);
sll_value = max(temp(temp < .9 * max(temp)));

sll = af_db_max - sll_value;
for i = hpbw_idx
    text(theta(i), af_db(i)+15, [sprintf("\\lceil%.2f\\pi", theta(i)/pi); sprintf("\\lfloor%.3g", round(af_db(i),4))], "HorizontalAlignment", "left", "Color", [1,.7,.7]);
end
title(sprintf("Side Lobe level = %.2f db, \\alpha = %g", sll, alpha))



for i = 1:length(hpbw_idx)/2
    drawbrace([theta(hpbw_idx(2*i)), af_db(hpbw_idx(2*i))-5], [theta(hpbw_idx(2*i-1)), af_db(hpbw_idx(1))-5], 0.01*diff(theta(hpbw_idx(2*i-1:2*i))), "Color", "red");
    text(mean(theta(hpbw_idx(2*i-1:2*i))), mean(af_db(hpbw_idx(2*i-1:2*i)))-20, sprintf("HPBW = %.2f\\pi", diff(theta(hpbw_idx(2*i-1:2*i)))/pi), "HorizontalAlignment", "center", "Color", [1,0,0]);
end

exportgraphics(gcf, 'results/af-db-with-keiser.pdf', 'Append', false);

