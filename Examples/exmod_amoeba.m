% FILE: amoeba.m
% trapezoid (see page 14 in: Tiago Anselmo et al 2020 J. Phys. A: Math.
% Theor. 53, 355201)
% 27-4-2025
clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n   =  2^12;
t   = (0:2*pi/n:2*pi-2*pi/n).';
%% 
s1    =  sin(t);
s2    =  sin(2.*t);
s4    =  sin(4.*t);
c1    =  cos(t);
c2    =  cos(2.*t);
es    =  exp(sin(t));
ec    =  exp(cos(t));
R     =  ec.*c2.^2+es.*s2.^2;
Rp    = -s1.*c2.^2.*ec+c1.*s2.^2.*es+2.*s4.*(es-ec);
et    =   R.*exp(i.*t);
etp   =  (Rp+i.*R).*exp(i.*t);
%
vz    = [et(1) ; et(n/4+1) ; et(2*n/4+1) ; et(3*n/4+1)];
%
alpha =  0;
[zet,zetp] = mapdisk(et,etp,n,alpha,'b');
vw = [zet(1) ; zet(n/4+1) ; zet(2*n/4+1) ; zet(3*n/4+1)];
modn = moddisk(vw(1),vw(2),vw(3),vw(4))
% 
%%
vs  = [0 ; 1   ; 1+modn*i  ; modn*i     ];
%
[ets,etsp]=polygonp(vs,n/4);
alphas = (1+modn*i)/2;
[zets,zetsp] = mapdisk(ets,etsp,n,alphas,'b');
vws = [zets(1) ; zets(n/4+1) ; zets(2*n/4+1) ; zets(3*n/4+1)];
%%
[x , y ] =  meshgrid(linspace(-1.5,2.8,1000),linspace(-1.2,2,1000));
zm       =  x+i*y;
zo       =  zm(:);
to       =  angle(zo);
ro       = (exp(cos(to)).*cos(2*to).^2+exp(sin(to)).*sin(2*to).^2);
zo(abs(zo)>=ro)  =  NaN+i*NaN;
%
z       =  zo(abs(zo)>=0).';
hw      =  fcau(et,etp,zet,z);
tw      =  mobius(hw,vw(1:3),vws(1:3));
w       =  fcau(zets,zetsp,ets,tw);
wo      =  NaN(size(zo));
wo(abs(zo)>=0) = w;
wm      =  NaN(size(zm));
wm(:) = wo;
%%
cntv  =  [0.00:0.05:0.95];
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
% axis square
axis([-1.5 2.8 -1.5 2.5])
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))
% axis square
% [cnt1,cnt2] = contour(real(zm),imag(zm),real(wm),cntv,'r','LineWidth',1.5);
contourf(real(zm),imag(zm),real(wm),cntv)
colormap jet
colr = 0.8.*colormap+0.2.*ones(size(colormap));
brighten(colr,0.4)
caxis([0  1])
colorbar 
crv = et(1:n);
plot(real(crv),imag(crv),'b','LineWidth',2)
plot(real(vz),imag(vz),'sk','MarkerSize',8,'MarkerFaceColor','k');
% plot(real(alpha),imag(alpha),'sk','MarkerSize',8,'MarkerFaceColor','k');
set(gca,'XTick',[-2:1:3]);
set(gca,'YTick',[-2:1:3]);
set(fig1,'PaperSize',[6  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_mdamoba.pdf', '-dpdf', '-fillpage');
%%
figure
hold on; box on
crv = zet(1:n);
% plot(real(hw),imag(hw),'.r');
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(0,0,'pr','MarkerSize',10,'MarkerFaceColor','r');
plot(real(vw),imag(vw),'dr','MarkerSize',8,'MarkerFaceColor','r');
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',18)
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))   
axis([-1.01 1.01 -1.01 1.01])
set(gca,'FontSize',18)
%%
figure
hold on; box on
crv = zets(1:n);
% plot(real(tw),imag(tw),'.r');
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(0,0,'pr','MarkerSize',10,'MarkerFaceColor','r');
plot(real(vws),imag(vws),'dr','MarkerSize',8,'MarkerFaceColor','r');
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',18)
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))   
axis([-1.01 1.01 -1.01 1.01])
set(gca,'FontSize',18)
%%
figure
hold on; box on
crv = ets(1:n);
plot(real(crv),imag(crv),'b','LineWidth',1.5)
plot(real(alphas),imag(alphas),'pk','MarkerSize',10,'MarkerFaceColor','k');
plot(real(vs),imag(vs),'sk','MarkerSize',8,'MarkerFaceColor','k');
grid on; grid('minor')
set(gca, 'XMinorTick','on'); set(gca, 'YMinorTick','on')
ax=gca; ax.GridAlpha=0.5; ax.MinorGridAlpha=0.5;
set(gca,'FontSize',18)
axis equal
set(gca,'LooseInset',get(gca,'TightInset'))   
axis([-0.2 1.2 -0.2 3.2])
set(gca,'FontSize',18)
%%