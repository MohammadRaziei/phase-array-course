KB = 1.38e-23;
T0 = 290;
BW = 10e6;
F_db = 5;
F = db2pow(F_db);

N0 = KB * T0 * BW * F
N0_db = pow2db(N0)

snr_min_db = 8;

p_min_db = snr_min_db + N0_db;
p_min = db2pow(p_min_db);
sigma = .5;
c = 3e8;
f = 8e9;
lambda = c / f;
P0 = 2;
G0_db = 6;
G0 = db2pow(G0_db);
% \frac{P_{\text{min}}(4 \pi)^3 R_{\text{max}}^4}{ P_0 G_0^2 \lambda^2 \sigma}
R_max = 100e3;
N3 = (p_min * (4*pi)^3 * R_max^4)/(P0 * G0^2 * lambda^2 * sigma);
N = ceil(N3^(1/3));

Gl = 7.61e-10;
G_max = G0 * 20

SLL = Gl / G_max
