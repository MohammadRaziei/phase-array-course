%% Q2
clear; close all; clc;
mkdir results
addpath ../common/

current_dir_part = strsplit(pwd, filesep()); % split path into parts
trogon_path = strjoin([current_dir_part(1:end-1), 'HW', 'Trogon.csv'], filesep()); % rebuild path excluding last part


data = readtable(trogon_path, "VariableNamingRule", "preserve");
disp(data.Properties.VariableNames)
data.Properties.VariableNames = {'Phi','Theta','dB'};

N = length(data.Variables);
Phi = reshape(data.Phi, sqrt(N), sqrt(N));
Theta = reshape(data.Theta, sqrt(N), sqrt(N));
DB = reshape(data.dB, sqrt(N), sqrt(N));

max(diff(Phi,1,1), [], "all"), min(Phi, [], "all"), max(Phi, [], "all")
max(diff(Theta,1,1), [], "all"), min(Theta, [], "all"), max(Theta, [], "all")


phiRad = deg2rad(Phi);
thetaRad = deg2rad(Theta);




figure;
surf(phiRad, thetaRad, DB)
shading interp
set(gca, "view", [-60, 55])
colorbar
xlabel("\phi (rad)")
ylabel("\theta (rad)")
zlabel("dB")
title('3D plot over angles');
exportgraphics(gcf, 'results/3d-plot-angles.pdf', 'Append', false);




figure;
imagesc(unique(phiRad), unique(thetaRad), DB)
xlabel("\phi (rad)")
ylabel("\theta (rad)")
title('2D plot over angles');
colorbar
exportgraphics(gcf, 'results/2d-plot-angles.pdf', 'Append', false);



[x,y,z] = sph2cart(phiRad, thetaRad, 10.^(DB/20));


figure
surf(x,y,z)
colormap([.5,.5,.5])
title('Spatial Antenna Pattern (non-dB scale)');
set(gca, "view", [-15, 15])
xlabel("x")
ylabel("y")
zlabel("z")
exportgraphics(gcf, 'results/spatial-antenna-pattern.pdf', 'Append', false);

%%
db_phi_0 = DB(Phi==0);
theta_phi_0 = Theta(Phi==0);
[db_max, db_max_idx] = max(db_phi_0);
db_hpbw = db_max - 3;
x = find((db_phi_0 - db_hpbw) >= 0);
hpbw_idx_3db = [min(x) max(x)+1];

figure
polarplot(deg2rad(theta_phi_0), 10.^(db_phi_0/20))
set(gca, 'ThetaZeroLocation', "top");
title('Antenna Pattern (where \phi = 0)')
exportgraphics(gcf, 'results/phi0-polar.pdf', 'Append', false);

figure('units','normalized','position',[0 .25 1 .5])

hold on
plot(theta_phi_0, db_phi_0)
xlim([-1, 1]*180);
ylim([-35, 20])
yline(db_hpbw);
plot(theta_phi_0(hpbw_idx_3db), db_phi_0(hpbw_idx_3db), 'ob')
ylabel('Power (dB)')
xlabel('\theta (degree)')
title('Power per \theta  (where \phi = 0)')

text(theta_phi_0(hpbw_idx_3db(1))-5, db_phi_0(hpbw_idx_3db(1))+2, sprintf("%.2g^\\circ", theta_phi_0(hpbw_idx_3db(1))), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
text(theta_phi_0(hpbw_idx_3db(2))+8, db_phi_0(hpbw_idx_3db(2))+2, sprintf("%.2g^\\circ", theta_phi_0(hpbw_idx_3db(2))), "HorizontalAlignment", "center", "Color", .3*ones(1,3))

drawbrace([theta_phi_0(hpbw_idx_3db(2)),db_phi_0(hpbw_idx_3db(2))-2], [theta_phi_0(hpbw_idx_3db(1)),db_phi_0(hpbw_idx_3db(1))-2], .005, "Color", "blue")
text(mean(theta_phi_0(hpbw_idx_3db)), mean(db_phi_0(hpbw_idx_3db))-5, sprintf("HPBW = %g^\\circ", diff(theta_phi_0(hpbw_idx_3db))), "HorizontalAlignment", "center", "Color", .3*ones(1,3))



[x,minIdx] = findpeaks(-db_phi_0);
nullIdx = minIdx(-x<0).';
[x,maxIdx] = findpeaks(db_phi_0);
maxIdx = [maxIdx(x<0).' db_max_idx];


plot(theta_phi_0(nullIdx), db_phi_0(nullIdx), 'or')
plot(theta_phi_0(maxIdx), db_phi_0(maxIdx), 'xr')

for i = nullIdx
    text(theta_phi_0(i)-5, db_phi_0(i)-2, sprintf("(%g\\circ, %.4g)", theta_phi_0(i), db_phi_0(i)), "HorizontalAlignment", "left", "Color", .3*ones(1,3))
end

for i = maxIdx
    text(theta_phi_0(i), db_phi_0(i)+2, sprintf("(%g\\circ, %.4g)", theta_phi_0(i), db_phi_0(i)), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
    if i == db_max_idx
        text(theta_phi_0(i), db_phi_0(i)+6, sprintf("antenna peak gain = %.4g dB", db_phi_0(i)), "HorizontalAlignment", "center", "Color", "blue")
    else
        text(theta_phi_0(i), db_phi_0(i)-8, sprintf("SLL = %.4g dB", db_max - db_phi_0(i)), "HorizontalAlignment", "center", "Color", "blue")
    end
end



exportgraphics(gcf, 'results/phi0.pdf', 'Append', false);

%%
db_phi_90 = DB(Phi==90);
theta_phi_90 = Theta(Phi==90);
[db_max, db_max_idx] = max(db_phi_90);
db_hpbw = db_max - 3;
x = find((db_phi_90 - db_hpbw) >= 0);
hpbw_idx_3db = [min(x) max(x)+1];

figure
polarplot(deg2rad(theta_phi_90), 10.^(db_phi_90/20))
set(gca, 'ThetaZeroLocation', "top");
title('Antenna Pattern (where \phi = 0)')
exportgraphics(gcf, 'results/phi90-polar.pdf', 'Append', false);

figure('units','normalized','position',[0 .25 1 .5])

hold on
plot(theta_phi_90, db_phi_90)
xlim([-1, 1]*180);
ylim([-35, 20])
yline(db_hpbw);
plot(theta_phi_90(hpbw_idx_3db), db_phi_90(hpbw_idx_3db), 'ob')
ylabel('Power (dB)')
xlabel('\theta (degree)')
title('Power per \theta  (where \phi = 90)')

text(theta_phi_90(hpbw_idx_3db(1))-5, db_phi_90(hpbw_idx_3db(1))+2, sprintf("%.2g^\\circ", theta_phi_90(hpbw_idx_3db(1))), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
text(theta_phi_90(hpbw_idx_3db(2))+8, db_phi_90(hpbw_idx_3db(2))+2, sprintf("%.2g^\\circ", theta_phi_90(hpbw_idx_3db(2))), "HorizontalAlignment", "center", "Color", .3*ones(1,3))

drawbrace([theta_phi_90(hpbw_idx_3db(2)),db_phi_90(hpbw_idx_3db(2))-2], [theta_phi_90(hpbw_idx_3db(1)),db_phi_90(hpbw_idx_3db(1))-2], .005, "Color", "blue")
text(mean(theta_phi_90(hpbw_idx_3db)), mean(db_phi_90(hpbw_idx_3db))-5, sprintf("HPBW = %g^\\circ", diff(theta_phi_90(hpbw_idx_3db))), "HorizontalAlignment", "center", "Color", .3*ones(1,3))



[x,minIdx] = findpeaks(-db_phi_90);
nullIdx = minIdx(-x<0).';
[x,maxIdx] = findpeaks(db_phi_90);
maxIdx = [maxIdx(x<0).' db_max_idx];


plot(theta_phi_90(nullIdx), db_phi_90(nullIdx), 'or')
plot(theta_phi_90(maxIdx), db_phi_90(maxIdx), 'xr')

for i = nullIdx
    text(theta_phi_90(i)-5, db_phi_90(i)-2, sprintf("(%g\\circ, %.4g)", theta_phi_90(i), db_phi_90(i)), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
end

for i = maxIdx
    text(theta_phi_90(i), db_phi_90(i)+2, sprintf("(%g\\circ, %.4g)", theta_phi_90(i), db_phi_90(i)), "HorizontalAlignment", "center", "Color", .3*ones(1,3))
    if i == db_max_idx
        text(theta_phi_90(i), db_phi_90(i)+6, sprintf("antenna peak gain = %.4g dB", db_phi_90(i)), "HorizontalAlignment", "center", "Color", "blue")
    else
        text(theta_phi_90(i), db_phi_90(i)-8, sprintf("SLL = %.4g dB", db_max - db_phi_90(i)), "HorizontalAlignment", "center", "Color", "blue")
    end
end



exportgraphics(gcf, 'results/phi90.pdf', 'Append', false);

%%

N = 4;
lambda = 3e8 / 60.48e9;
d_lambda_x = 7e-3 / lambda;
d_lambda_y = 5e-3 / lambda;
AF_x =@(theta, phi) sin((N-1) * pi * d_lambda_y * sin(theta) .* cos(phi) )./ sin( pi * d_lambda_y * sin(theta) .* cos(phi) ); 
AF_y =@(theta, phi) sin((N-1) * pi * d_lambda_y * sin(theta) .* sin(phi) )./ sin( pi * d_lambda_y * sin(theta) .* sin(phi) ); 


AF_xy = AF_x(thetaRad, phiRad) .* AF_y(thetaRad, phiRad);
AF_total = DB + mag2db(AF_xy);


figure;
imagesc(unique(phiRad), unique(thetaRad), mag2db(AF_xy))
xlabel("\phi (rad)")
ylabel("\theta (rad)")
title('2D plot over angles');
colorbar
exportgraphics(gcf, 'results/2d-plot-angles-4x4.pdf', 'Append', false);


[x,y,z] = sph2cart(phiRad, thetaRad, AF_xy);
figure
surf(x,y,z)
colormap([.5,.5,.5])
title('Spatial Antenna Pattern (non-dB scale)');
set(gca, "view", [-15, 15])
xlabel("x")
ylabel("y")
zlabel("z")
exportgraphics(gcf, 'results/spatial-antenna-pattern-4x4.pdf', 'Append', false);


figure;
surf(phiRad, thetaRad, mag2db(AF_xy))
shading interp
set(gca, "view", [-60, 55])
colorbar
xlabel("\phi (rad)")
ylabel("\theta (rad)")
zlabel("dB")
title('3D plot over angles');
exportgraphics(gcf, 'results/3d-plot-angles-4x4.pdf', 'Append', false);



%%


figure;
imagesc(unique(phiRad), unique(thetaRad), AF_total)
xlabel("\phi (rad)")
ylabel("\theta (rad)")
title('2D plot over angles');
colorbar
exportgraphics(gcf, 'results/2d-plot-angles-4x4-total.pdf', 'Append', false);




figure;
surf(phiRad, thetaRad, AF_total)
shading interp
set(gca, "view", [-60, 55])
colorbar
xlabel("\phi (rad)")
ylabel("\theta (rad)")
zlabel("dB")
title('3D plot over angles');
exportgraphics(gcf, 'results/3d-plot-angles-4x4-total.pdf', 'Append', false);



[x,y,z] = sph2cart(phiRad, thetaRad, 10.^(AF_total/20));


figure
surf(x,y,z)
colormap([.5,.5,.5])
title('Spatial Antenna Pattern (non-dB scale)');
set(gca, "view", [-15, 15])
xlabel("x")
ylabel("y")
zlabel("z")
exportgraphics(gcf, 'results/spatial-antenna-pattern-4x4-total.pdf', 'Append', false);