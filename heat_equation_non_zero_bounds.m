clc;
clear;
close all;

%aluminium
density=2700;%kg/m3
heat_cap=897;%J/(kg*K)
thermal_cond=237;%W/(m·K)
k=thermal_cond/(heat_cap*density);%m^2/c

N = 1000;%stała, do której liczymy końcowe równanie
L=10;%długość walca
dx = 0.1;%krok z jakim liczymy wartości dla odległości
dt = 0.1;%krok z jakim liczymy wartości dla czasu
X = 0:dx:L;%wektor odległości
T = (0:dt:L)';%wektor czasu
u = zeros(length(T),length(X));%szukana funkcja temperatury
fx=@(x) x.*(x.^2-3.*L.*x+2.*L.^2);%warunek początkowy
T0=0;
TL=100;
uEx=@(x) T0+(TL-T0)*x/TL;%"equilibrium temperature"

for n=1:1:N
    expr=@(x) (fx(x)-uEx(x)).*sin(n.*pi.*x./L);
    Bn=(2/L)*integral(fx,0,L);
    u=u+Bn.*sin(n.*pi.*X./L).*exp(-k.*(n.*pi./L).^2.*T);
end

% Rysujemy wykres 3D
surf(X,T,u)
xlabel('distance','Fontsize',20);
ylabel('time','Fontsize',20);
zlabel('temperature u(x,t)','Fontsize',20);
