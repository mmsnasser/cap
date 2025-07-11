% FILE: trapezoid.m
% trapezoid (see page 14 in: Tiago Anselmo et al 2020 J. Phys. A: Math.
% Theor. 53, 355201)
% 24-4-2025
clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n   =  2^12;
t   = (0:2*pi/n:2*pi-2*pi/n).';
% 
%%
[x , y ] =  meshgrid(0.1:0.05:3);
zm       =  x+i*y;
zm(x+y<=1) = NaN+i*NaN;
zo       =  zm(:);
z       =  zo(abs(zo)>=0).';
%%
for kk=1:length(z)
v1  =  0; v2  =  1; v3  =  z(kk); v4  = i ;
%
vz  = [v1  ; v2  ; v3  ; v4   ];
%
% Find the exact value of the modulus
mode(kk,1) =  QM(v3,v4);
%
% approximate modulus
[et,etp]=polygonp(vz,n/4);
%
alpha = (v2+v4)/2;
[zet,zetp] = mapdisk(et,etp,n,alpha,'b');
vw = [zet(1) ; zet(n/4+1) ; zet(2*n/4+1) ; zet(3*n/4+1)];
modn(kk,1) = moddisk(vw(1),vw(2),vw(3),vw(4));
% 
rerr(kk,1) = abs(mode(kk)-modn(kk))/mode(kk);
% 
end
%%
modno      =  NaN(size(zo));
modno(abs(zo)>=0) = modn;
modnv      =  NaN(size(zm));
modnv(:)   = modno;
%
rerro      =  NaN(size(zo));
rerro(abs(zo)>=0) = rerr;
rerrv      =  NaN(size(zm));
rerrv(:)   = rerro;
%%
cntv  =  [0.1:0.1:3];
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
% [cnt1,cnt2] = contour(real(zm),imag(zm),real(wm),cntv,'r','LineWidth',1.5);
contourf(real(zm),imag(zm),modnv,cntv)
colormap jet
colr = 0.8.*colormap+0.2.*ones(size(colormap));
brighten(colr,0.4)
caxis([0  3])
colorbar 
axis([0 3 0 3])
set(gca,'XTick',[0:1:3]);
set(gca,'YTick',[0:1:3]);
set(fig1,'PaperSize',[6  5]);
grid on; axis square
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
%
print(fig1, 'fig_mdconv.pdf', '-dpdf', '-fillpage');
%%
fig2=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
surf(real(zm),imag(zm),rerrv)
% zscale log
colormap jet
shading interp
caxis([0 1e-10])
box on
set(fig2,'PaperSize',[6  5]);
% axis([0 3 0 3 0 1.05])
view([60 30]) 
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-1:1:8],'FontSize',18);
set(gca,'YTick',[0:1:6]);
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gcf,'Renderer','zbuffer')
print(fig2, 'fig_mdrerr.pdf', '-dpdf', '-fillpage');
%%