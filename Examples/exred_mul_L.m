clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n     =  2^11;
t     = (0:2*pi/n:2*pi-2*pi/n).';
%
cent  =  [0     ; 0.6  ; -0.3+0.5i ; -0.3-0.5i];
rad   =  [0.15  ; 0.2  ;  0.2      ;  0.2     ];
%%
m  =  length(rad);
et = [];  etp = [];
et    =    exp(i*t);
etp   =  i*exp(i*t);
for k=1:m 
    et  = [et  ; cent(k)+rad(k)*exp(-i*t)];
    etp = [etp ;      -i*rad(k)*exp(-i*t)];
end
%
%%
x = linspace(-1,1,200); 
y = linspace(-1,1,200); 
[xx,yy] = meshgrid(x,y); 
zz = xx+1i*yy;
zv = zz(:);
zv(abs(zv)>1-0.005) = NaN+i*NaN;
for k=1:m
    zv(abs(zv-cent(k))<rad(k)+0.005) = NaN+i*NaN;
end
z  = zv(abs(zv)>=0); z=z(:).';
%%
for k=1:length(z)
    k
    alpha   =  z(k);
    [~,~,c] =  diskcirslitmap(et,etp,n,alpha);
    redvp(k) = -log(c)/(2*pi);
end
% 
redv = NaN(size(zv));
redv(abs(zv)>=0)=redvp;
redm = reshape(redv,size(zz));
%%
cntv  =  [-0.35,-0.26,-0.2,-0.15,-0.12,-0.1,-0.09,-0.086];
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
contour(xx,yy,redm,cntv,'color','k','linewidth',1);
% plot(real(z),imag(z),'.r')
for k=1:m+1
    Jk = (k-1)*n+1:k*n;
    plot(real(et(Jk)),imag(et(Jk)),'b','LineWidth',2);
end
axis equal
axis([-1  1 -1  1])
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick',[-3:1:3]);
% set(gca,'YTick',[-3:1:3]);
set(fig1,'PaperSize',[5  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_mulL.pdf', '-dpdf', '-fillpage');
%%

