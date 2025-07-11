clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n     =  3*2^14;  
t     = (0:2*pi/n:2*pi-2*pi/n).';
%
ver      = [2i;3i;-2+3i;-2-2i;3-2i;3+i;i;0;2;2-i;-1-i;-1+2i];
[et,etp] =  polygonp(ver,n/length(ver));
%
alpha =  0.5-1.5i;
[zet,zetp,~] = mapdisk(et,etp,n,alpha,'b');
%
rho    = @(a,b)2*asinh(abs(a-b)./sqrt((1-abs(a).^2).*(1-abs(b).^2)));
inpoly = @(z,w) inpolygon(real(z),imag(z),real(w),imag(w));
%%
figure(1);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(real(zet),imag(zet),'k','LineWidth',2)
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-3:1:3]);
set(gca,'YTick',[-3:1:3]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
figure(2);
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(real(et),imag(et),'k','LineWidth',2)
axis([-2 3 -2 3])
axis square
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[-3:1:3]);
set(gca,'YTick',[-3:1:3]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
%
tvp = [0.01:0.01:0.99,0.991:0.001:0.999,0.9991:0.0001:0.9999,...
      0.99991:0.00001:0.99999,0.999991:0.000001:0.999999,...
      0.9999991:0.0000001:0.9999999,0.99999991:0.00000001:0.99999999,...
      0.999999991:0.000000001:0.999999999,0.9999999991:0.0000000001:0.9999999999,...
      0.99999999991:0.00000000001:0.99999999999,0.999999999991:0.000000000001:0.999999999999];
%
tv = [-tvp(end:-1:1),0,tvp];
for k=1:length(ver)
    a   = zet(((2*k-1)*n/length(ver))/2+1);
    La  = tv*a;
    Laz = fcau(zet,zetp,et,La);
    %
    figure(1)
    plot(real(La),imag(La),'LineWidth',1.5)
    %
    figure(2)
    plot(real(Laz),imag(Laz),'LineWidth',1.5)
end
%%
fig1=figure(2);
set(fig1,'PaperSize',[5  5]);
print(fig1, 'fig_hypgeod1.pdf', '-dpdf', '-fillpage');
%%
zp = [-1+2.5i,-1.5,0.5-1.5i,2.5-0.5i,1+0.5i];
wp = fcau(et,etp,zet,zp);
%%
figure(2)
plot(real(zp),imag(zp),'or','MarkerFaceColor','r');
%%
for k=1:length(zp)
    for j=1:length(zp)
        if k==j
            hdis(k,j) = 0;
        else
            hdis(k,j) = rho(wp(k),wp(j));
        end
    end
end
hdis
%%