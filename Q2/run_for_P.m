function run_for_P(P, prefix)

lambda = 1;
d = lambda/2;

assert((min(pdist(P))-d)<1e-5)

AF = @(P, lambda, theta, phi) abs(cell2mat(arrayfun(@(phi) arrayfun(@(theta) sum(exp(1j*2*pi/lambda *[sin(theta)* cos(phi), sin(theta)* sin(phi), cos(theta)]*P')), theta,'UniformOutput',true)', phi,'UniformOutput',false)));

phi = 0:0.004:2*pi;
theta = -pi:0.004:pi;


[Phi, Theta]=meshgrid(phi, theta);


%% a

af = AF(P, lambda, theta, phi);


figure;
plot3(P(:, 1), P(:, 2), P(:, 3), 'o', 'MarkerSize', 6, 'MarkerFaceColor', 'b');
grid on;
exportgraphics(gcf, sprintf('results/%s-p-plot3.pdf', prefix), 'Append', false);


figure;
imagesc(phi, theta, af)
xlabel("\phi"); ylabel("\theta");
exportgraphics(gcf, sprintf('results/%s-imagesc.pdf', prefix), 'Append', false);


figure
[x, y, z] = sph2cart(Phi, pi/2-Theta, af);
surf(x, y, z,  'EdgeColor', 'black', "EdgeAlpha", .05); colormap(.7*ones(1,3));
xlabel("x"); ylabel("y"); zlabel("z")
lim = [-30,30]; xlim(lim); ylim(lim); zlim(lim);

exportgraphics(gcf, sprintf('results/%s-cart.pdf', prefix), 'Append', false);



theta_phi_0 = Theta(Phi==0);
af_phi_0 = af(Phi==0);


[~,MaxIdx] = findpeaks(af_phi_0);
[~,MinIdx] = findpeaks(-af_phi_0);




figure('units','normalized','position',[0 .25 1 .5]); hold on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];



plot(theta_phi_0, af_phi_0)
xlabel("\theta"); ylabel("Magnitude")
xlim([-1, 1]*pi); xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))
ylim([-5, 30]);


plot(theta_phi_0(MaxIdx), af_phi_0(MaxIdx), ".b");
plot(theta_phi_0(MinIdx), af_phi_0(MinIdx), ".b");


for i = MaxIdx'
    text(theta_phi_0(i), af_phi_0(i)+1.8, sprintf("\\lceil%.2g\\pi", round(theta_phi_0(i)/pi, 2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
    text(theta_phi_0(i), af_phi_0(i)+.8, sprintf("\\lfloor%.2g", round(af_phi_0(i),2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
end
for i = MinIdx'
    text(theta_phi_0(i), af_phi_0(i)-.8, sprintf("\\lceil%.2g\\pi", round(theta_phi_0(i)/pi, 2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
    text(theta_phi_0(i), af_phi_0(i)-1.8, sprintf("\\lfloor%.2g", round(af_phi_0(i),2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
end


exportgraphics(gcf, sprintf('results/%s-phi-0-mag.pdf', prefix), 'Append', false);




figure('units','normalized','position',[0 .25 1 .5]); hold on
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3);
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];


af_phi_0_db = mag2db(af_phi_0);

plot(theta_phi_0, af_phi_0_db)

% ll = length(af_phi_0_db);
% is = floor(ll/4); ie = ceil(ll*3/4);
b = find(af_phi_0_db >= (max(af_phi_0_db) - 3));
a = diff([0; b]);
a = find((circshift(a,-1)>1) + (a>1));
if b(1) == 1, a = a(2:end-1); end 
idx_hpbw = reshape(b(a), 2, []); % 2 * n; each col each pair; 

xlim([-1, 1]*pi); xticks((-1:0.25:1)*pi)
xticklabels(arrayfun(@(s) sprintf("%g \\pi", s), (-1:0.25:1),'UniformOutput',true))
xlabel("\theta"); ylabel("dB")
ylim([-50, 40]);

plot(theta_phi_0(idx_hpbw(:)), af_phi_0_db(idx_hpbw(:)), "or");
yline((max(af_phi_0_db) - 3))

for i = 1:size(idx_hpbw, 2)
    text(theta_phi_0(idx_hpbw(1,i))-.2, af_phi_0_db(idx_hpbw(1,i))+2.5+2.5, sprintf("\\lceil%.3g\\pi", round(theta_phi_0(idx_hpbw(1,i))/pi, 3)), "HorizontalAlignment", "left", "Color", "red")
    text(theta_phi_0(idx_hpbw(1,i))-.2, af_phi_0_db(idx_hpbw(1,i))+2.5, sprintf("\\lfloor%.3g", round(af_phi_0_db(idx_hpbw(1,i))/pi, 3)), "HorizontalAlignment", "left", "Color", "red")
    text(theta_phi_0(idx_hpbw(2,i))+.02, af_phi_0_db(idx_hpbw(2,i))+2.5+2.5, sprintf("\\lceil%.3g\\pi", round(theta_phi_0(idx_hpbw(2,i))/pi, 3)), "HorizontalAlignment", "left", "Color", "red")
    text(theta_phi_0(idx_hpbw(2,i))+.02, af_phi_0_db(idx_hpbw(2,i))+2.5, sprintf("\\lfloor%.3g", round(af_phi_0_db(idx_hpbw(2,i))/pi, 3)), "HorizontalAlignment", "left", "Color", "red")

    text(mean(theta_phi_0(idx_hpbw(:,i))), mean(af_phi_0_db(idx_hpbw(:,i)))+12, sprintf("HPBW = %.3f\\pi", diff(theta_phi_0(idx_hpbw(:,i))/pi)), "HorizontalAlignment", "center", "Color", "blue")
end

plot(theta_phi_0(MaxIdx), af_phi_0_db(MaxIdx), ".b");
plot(theta_phi_0(MinIdx), af_phi_0_db(MinIdx), ".b");

for i = MaxIdx'
    text(theta_phi_0(i), af_phi_0_db(i)+2.5+2.5, sprintf("\\lceil%.2g\\pi", round(theta_phi_0(i)/pi, 2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
    text(theta_phi_0(i), af_phi_0_db(i)+2.5, sprintf("\\lfloor%.2g", round(af_phi_0_db(i),2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
end
for i = MinIdx'
    text(theta_phi_0(i), af_phi_0_db(i)-2.5, sprintf("\\lceil%.2g\\pi", round(theta_phi_0(i)/pi, 2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
    text(theta_phi_0(i), af_phi_0_db(i)-2.5-2.5, sprintf("\\lfloor%.2g", round(af_phi_0_db(i),2)), "HorizontalAlignment", "left", "Color", .4*ones(1,3))
end



exportgraphics(gcf, sprintf('results/%s-phi-0-db.pdf', prefix), 'Append', false);




