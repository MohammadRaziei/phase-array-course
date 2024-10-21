clear; close all; clc;
mkdir results

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


[x,y,z] = sph2cart(phiRad, thetaRad, 10.^(DB/20));

surf(x,y,z)
colormap([.5,.5,.5])


% colormap(jet);
% colorbar;
title('Surface Plot');

