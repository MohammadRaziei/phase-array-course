function arrayfactor = Array_beam_cal(R,phases,Tnew,Pnew)

f0=60e9;
wavelength=3e8/f0;
k0=2*pi/wavelength;
d=wavelength/2;
dx=wavelength/2;
dy=dx;
Noise=0.01;
T=(Tnew)*pi/180; %T=(Tnew-90)*pi/180;
P=Pnew*pi/180;
r=[sin(T)*cos(P) sin(T)*sin(P) cos(T)];
AF=0;
Nel=length(R);

for in=1:Nel;
    wpi=phases(in);
%     pr=(compare(V(in)+eps,Vref))
%     dvi=V(in)-Vref(pr);
%    wai=(Am(pr)+dvi*(Am(pr+1)-Am(pr)))/max(Am);
     AF=AF+exp(-j*wpi)*exp(j*k0*sum(R(in,:).*r))+Noise*randn;
end
arrayfactor= abs(AF);   

% R2=[-1.5 0 0
% -1.05 0 0
% 0.05 0 0
% 0.95 0 0
% 2.2 0 0];
