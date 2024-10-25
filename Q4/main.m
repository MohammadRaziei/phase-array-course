%% Q4
clear; close all; clc;
mkdir results
addpath ../common/
addpath ../HW/

rng(1)

%%
R1=(2.5e-3)*...
    [3 0 0
    2 0 0
    1 0 0
    0 0 0
    -1 0 0
    -2 0 0
    -3 0 0];

[~, fig1, fig2, ~, ~] = run_array_beam(R1, 4);
exportgraphics(fig1, 'results/array-beam-polar-R1.pdf', 'Append', false);
exportgraphics(fig2, 'results/array-beam-R1.pdf', 'Append', false);

csvwrite('results/R1.csv', R1);

%%

R2=(5e-3)*...
    [3 0 0
    2 0 0
    1 0 0
    0 0 0
    -1 0 0
    -2 0 0
    -3 0 0];

[~, fig1, fig2, ~, ~] = run_array_beam(R2, 4);
exportgraphics(fig1, 'results/array-beam-polar-R2.pdf', 'Append', false);
exportgraphics(fig2, 'results/array-beam-R2.pdf', 'Append', false);

csvwrite('results/R2.csv', R2);


%%

R3=(5e-3)*...
    [3+0.5*randn() 0 0
    2+0.5*randn() 0 0
    1+0.5*randn() 0 0
    0+0.5*randn() 0 0
    -1+0.5*+randn() 0 0
    -2+0.5*randn() 0 0
    -3+0.5*randn() 0 0];


[~, fig1, fig2, ~, ~] = run_array_beam(R3, 2);
exportgraphics(fig1, 'results/array-beam-polar-R3.pdf', 'Append', false);
exportgraphics(fig2, 'results/array-beam-R3.pdf', 'Append', false);

csvwrite('results/R3.csv', R3);

%%
R3=(5e-3)*...
    [3+0.5*randn() 0 0
    2+0.5*randn() 0 0
    1+0.5*randn() 0 0
    0+0.5*randn() 0 0
    -1+0.5*+randn() 0 0
    -2+0.5*randn() 0 0
    -3+0.5*randn() 0 0];

csvwrite('results/R3-2.csv', R3);


[~, fig1, fig2, ~, ~] = run_array_beam(R3,2);
exportgraphics(fig1, 'results/array-beam-polar-R3-2.pdf', 'Append', false);
exportgraphics(fig2, 'results/array-beam-R3-2.pdf', 'Append', false);




%%

rng(3)
num_iter = 800;


peak_ratio_max = 0;
peak_max_min= inf;
for i = 1:num_iter
    R=(5e-3)*...
        [3+0.5*randn() 0 0
        2+0.5*randn() 0 0
        1+0.5*randn() 0 0
        0+0.5*randn() 0 0
        -1+0.5*+randn() 0 0
        -2+0.5*randn() 0 0
        -3+0.5*randn() 0 0];
    [arrayfactor, ~, ~, ~, MaxIdx] = run_array_beam(R,2, false);
    peaks = arrayfactor(MaxIdx);
    peak_max= max(peaks);
    peak_max2 = max(peaks(peaks<.95*peak_max));
    peak_ratio = peak_max / peak_max2;
    if peak_ratio_max < peak_ratio
        peak_ratio_max = peak_ratio;
        Ropt = R;
    end
    if peak_max_min > peak_max2
        peak_max_min = peak_max2;
        Ropt2 = R;
    end
end

disp(peak_ratio_max);

[arrayfactor, fig1, fig2, ~, MaxIdx] = run_array_beam(Ropt,2, true);
peaks = arrayfactor(MaxIdx);
peak_max= max(peaks);
peak_max2 = max(peaks(peaks<.95*peak_max));
peak_ratio = peak_max / peak_max2;

csvwrite('results/peak-opt.csv', [peak_max, peak_max2, peak_ratio]);
csvwrite('results/Ropt.csv', Ropt);
csvwrite('results/Ropt-R2.csv', Ropt-R2);

exportgraphics(fig1, 'results/array-beam-polar-Ropt.pdf', 'Append', false);
exportgraphics(fig2, 'results/array-beam-Ropt.pdf', 'Append', false);

disp(peak_max_min);

[arrayfactor, fig1, fig2, ~, MaxIdx] = run_array_beam(Ropt2,2, true);
peaks = arrayfactor(MaxIdx);
peak_max= max(peaks);
peak_max2 = max(peaks(peaks<.95*peak_max));
peak_ratio = peak_max / peak_max2;

csvwrite('results/peak-opt-2.csv', [peak_max, peak_max2, peak_ratio]);
csvwrite('results/Ropt-2.csv', Ropt2);
csvwrite('results/Ropt-R2-2.csv', Ropt2-R2);

exportgraphics(fig1, 'results/array-beam-polar-Ropt-2.pdf', 'Append', false);
exportgraphics(fig2, 'results/array-beam-Ropt-2.pdf', 'Append', false);



%%
