clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
tau   =  0.5;
% 
n     =  2^12;
t     = (0:2*pi/n:2*pi-2*pi/n).';
et    =    cosh(tau+i*t);
etp   =  i*sinh(tau+i*t);
%
%%
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
%
x = linspace(-cosh(tau),cosh(tau),100); 
y = linspace(-sinh(tau),sinh(tau),100); 
[xx,yy] = meshgrid(x,y); 
zz = xx+1i*yy;
zv = zz(:);
ii = inpoly(1.015*zv,et); 
zv(~ii)=NaN+i*NaN;
z  = zv(abs(zv)>=0); z=z(:).';
%%
for k=1:length(z)
    k
    alpha   =  z(k);
    [~,~,c] =  mapdisk(et,etp,n,alpha,'b');
    redvp(k) = -log(c)/(2*pi);
end
%% 
redv = NaN(size(zv));
redv(abs(zv)>=0)=redvp;
redm = reshape(redv,size(zz));
%%
cntv  =  [-0.4,-0.28,-0.2,-0.15,-0.118,-0.095,-0.082,-0.076,-0.0725];
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
contour(xx,yy,redm,cntv,'color','k','linewidth',1);
plot(real(et),imag(et),'b','LineWidth',2)
axis equal
axis([-cosh(tau)  cosh(tau) -sinh(tau)  sinh(tau)])
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick',[-3:1:3]);
% set(gca,'YTick',[-3:1:3]);
set(fig1,'PaperSize',[6  3.0]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_inellipseL.pdf', '-dpdf', '-fillpage');
%%

