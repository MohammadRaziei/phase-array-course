%% Q3
clear; close all; clc;
mkdir results
addpath ../common/
addpath ../HW/


%%


current_dir_part = strsplit(pwd, filesep()); % split path into parts
trogon_path = strjoin([current_dir_part(1:end-1), 'HW', 'Trogon_Sparam.csv'], filesep()); % rebuild path excluding last part


data = readtable(trogon_path, "VariableNamingRule", "preserve");
data.Properties.VariableNames = {'freq','s11','s12', 's21', 's22'};



%%
disp([min(data.freq) max(data.freq)])
figure
plot(data.freq, data.s11)
ylim([-35,0])
xlabel("Frequency (GHz)"); ylabel("S11 parameter (dB)")
exportgraphics(gcf, 'results/s11.pdf', 'Append', false);
figure
plot(data.freq, data.s12)
ylim([-35,0])
xlabel("Frequency (GHz)"); ylabel("S12 parameter (dB)")
exportgraphics(gcf, 'results/s12.pdf', 'Append', false);
figure
plot(data.freq, data.s21)
ylim([-35,0])
xlabel("Frequency (GHz)"); ylabel("S21 parameter (dB)")
exportgraphics(gcf, 'results/s21.pdf', 'Append', false);
figure
plot(data.freq, data.s22)
ylim([-35,0])
xlabel("Frequency (GHz)"); ylabel("S22 parameter (dB)")
exportgraphics(gcf, 'results/s22.pdf', 'Append', false);
%%
figure; hold on
plot(data.freq, data.s11, "LineWidth",1)
plot(data.freq, data.s12,"LineWidth",1)
plot(data.freq, data.s21, "LineWidth",1)
plot(data.freq, data.s22, "LineWidth",1)
legend("S_{11}", "S_{12}", "S_{21}", "S_{22}", 'Location', "southwest")
ylim([-35,0])
exportgraphics(gcf, 'results/S-param.pdf', 'Append', false);

%%
figure; hold on
plot(data.freq, data.s11)
ylim([-35,0])
yline(-10)

x = find((data.s11 +10) <= 0);
bw_idx = [min(x) max(x)];

plot(data.freq(bw_idx), data.s11(bw_idx), "ob")

for i = bw_idx
    text(data.freq(i), data.s11(i)+3, sprintf("%.3g", data.freq(i)), "HorizontalAlignment", "center", "Color", .4*ones(1,3))
end

drawbrace([data.freq(bw_idx(1)),data.s11(bw_idx(1))+5], [data.freq(bw_idx(2)),data.s11(bw_idx(2))+5], 5, "Color", "blue")
text(mean(data.freq(bw_idx)), mean(data.s11(bw_idx))+8, sprintf("BW = %g (Hz)", diff(data.freq(bw_idx))), "HorizontalAlignment", "center", "Color", .3*ones(1,3))

exportgraphics(gcf, 'results/s11-10db.pdf', 'Append', false);
