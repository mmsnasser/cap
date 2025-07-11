clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n     =  2^12;  
t     = (0:2*pi/n:2*pi-2*pi/n).';
et    = (2+cos(5*t)).*exp(i*t);
etp   = (-5*sin(5*t)+i*(2+cos(5*t))).*exp(i*t);
%
alpha =  1;
[zet,zetp,~] = mapdisk(et,etp,n,alpha,'b');
%
rho    = @(a,b)2*asinh(abs(a-b)./sqrt((1-abs(a).^2).*(1-abs(b).^2)));
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
%%
x = linspace(-3,3,200); 
y = linspace(-3,3,200);
[xx,yy] = meshgrid(x,y); 
zz = xx+1i*yy;
ii = inpoly(zz,et); 
zz(~ii)=NaN+i*NaN;
zv = zz(:);
z  = zv(abs(zv)>=0); z=z(:).';
% 
w = fcau(et,etp,zet,z);
wv = (1+i)*NaN(size(zv));
wv(abs(zv)>=0)=w;
ww = reshape(wv,size(zz));
uu = rho(ww,0);
%
%%
cntv  =  [0,0.5:1:9.5];
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
contourf(xx,yy,uu,cntv)
colormap jet
colr = 0.8.*colormap+0.2.*ones(size(colormap));
brighten(colr,0.4)
caxis([0  10])
% colorbar 
plot(real(et),imag(et),'b','LineWidth',2)
axis square
axis([-3 3 -3 3])
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-3:1:3]);
set(gca,'YTick',[-3:1:3]);
set(fig1,'PaperSize',[6  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_hypdisk2.pdf', '-dpdf', '-fillpage');
%%
