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

[arrayfactor, fig1, fig2, MinIdx, MaxIdx] = run_array_beam(R1, 4);
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

[arrayfactor, fig1, fig2, MinIdx, MaxIdx] = run_array_beam(R2, 4);
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


[arrayfactor, fig1, fig2, MinIdx, MaxIdx] = run_array_beam(R3, 4);
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


[arrayfactor, fig1, fig2, MinIdx, MaxIdx] = run_array_beam(R3,4);
exportgraphics(fig1, 'results/array-beam-polar-R3-2.pdf', 'Append', false);
exportgraphics(fig2, 'results/array-beam-R3-2.pdf', 'Append', false);









