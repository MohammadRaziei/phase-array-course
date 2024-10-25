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
plot(data.freq, data.s11, "LineWidth", 1)
ylim([-35,0])
xlabel("Frequency (GHz)"); ylabel("S11 parameter (dB)")
exportgraphics(gcf, 'results/s11.pdf', 'Append', false);
figure
plot(data.freq, data.s12, "LineWidth", 1)
ylim([-35,0])
xlabel("Frequency (GHz)"); ylabel("S12 parameter (dB)")
exportgraphics(gcf, 'results/s12.pdf', 'Append', false);
figure
plot(data.freq, data.s21, "LineWidth", 1)
ylim([-35,0])
xlabel("Frequency (GHz)"); ylabel("S21 parameter (dB)")
exportgraphics(gcf, 'results/s21.pdf', 'Append', false);
figure
plot(data.freq, data.s22, "LineWidth", 1)
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
xlabel("Frequency (GHz)"); ylabel("S-parameters (dB)")
exportgraphics(gcf, 'results/S-param.pdf', 'Append', false);

%%
[~, MinIdx] = findpeaks(-data.s11);

figure; hold on
plot(data.freq, data.s11, "LineWidth", 1)
ylim([-35,0])
yline(-10)
xlabel("Frequency (GHz)"); ylabel("S11 parameter (dB)")

x = find((data.s11 +10) <= 0);
bw_idx = [min(x) max(x)];

plot(data.freq(bw_idx), data.s11(bw_idx), "ob")
plot(data.freq(MinIdx), data.s11(MinIdx), "or")

f0 = data.freq(MinIdx);

for i = bw_idx
    text(data.freq(i), data.s11(i)+3, sprintf("%.3g", data.freq(i)), "HorizontalAlignment", "center", "Color", .4*ones(1,3))
end
for i = MinIdx
    text(data.freq(i), data.s11(i)-2, sprintf("(%.3g, %.3g)", data.freq(i), data.s11(i)), "HorizontalAlignment", "center", "Color", .4*ones(1,3))
    text(data.freq(i)+4, data.s11(i), sprintf("\\rightarrow f_0 = %.3g (GHz)", data.freq(i)), "HorizontalAlignment", "center", "Color", "blue")
end



drawbrace([data.freq(bw_idx(1)),data.s11(bw_idx(1))+5], [data.freq(bw_idx(2)),data.s11(bw_idx(2))+5], 5, "Color", "blue")
text(mean(data.freq(bw_idx)), mean(data.s11(bw_idx))+8, sprintf("BW = %g (GHz) ,  BW_{frac} = %.4g", diff(data.freq(bw_idx)), diff(data.freq(bw_idx))/f0), "HorizontalAlignment", "center", "Color", "blue")

exportgraphics(gcf, 'results/s11-10db.pdf', 'Append', false);
%%
figure; hold on
plot(data.freq, data.s12, "LineWidth", 1)
ylim([-35,0])
xline(f0)
plot(data.freq(MinIdx), data.s12(MinIdx), "or")
xlabel("Frequency (GHz)"); ylabel("S12 parameter (dB)")

text(data.freq(MinIdx)-.5, data.s12(MinIdx)-2, sprintf("Antenna coupling at the f_0: %.3g", data.s12(i)), "HorizontalAlignment", "right", "Color", .3*ones(1,3))
text(data.freq(MinIdx)-.5, data.s12(MinIdx)-5, sprintf("where f_0= %.3g", data.freq(i)), "HorizontalAlignment", "right", "Color", .3*ones(1,3))
exportgraphics(gcf, 'results/s12-f0.pdf', 'Append', false);
