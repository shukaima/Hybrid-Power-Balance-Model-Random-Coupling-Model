%%
clc
marker = 777;
nom = 140;
% Gr=real(Ydiag(1:nom,1:nom,marker));
% Br=imag(Ydiag(1:nom,1:nom,marker));
Gr=grad(1:nom,1:nom,marker);
Br=brad(1:nom,1:nom,marker);
Yd=Gr+1i*Br;
Zd=inv(Yd);
ZdH=0.5*(Zd+Zd');
YdH=0.5*(Yd+Yd');
aa=(Yd.')';
% ZdH=Zd';
% YdH=Yd';
qq=abs(trace(ZdH*Yd*ZdH*(Yd.')'));
freq=75*10^9+marker/10001*(110-75)*10^9;
c=3*10^8;
k=2*pi*freq/c;

sig1 = 2*pi/k^2*qq

sig2 = pi*(6.35*10^-3)^2/2
%%
Acav=0.016;
freq=100*10^9;
c=3*10^8;
lambda=c/freq;
mu0=4*pi*10^-7;
sigma=6*10^7;
sigwall=4*pi*Acav/3/lambda*sqrt(1/pi/freq/mu0/sigma);
%%
Pwall=0.068;
%Pwall=0.0193;
Sc=Pwall/sigwall;
sig3=0.093/Sc
%sig3=0.027/Sc
%%
clc
sig1
sig2
sig3