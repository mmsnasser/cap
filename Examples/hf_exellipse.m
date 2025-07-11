clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
tau   =  0.9;
zo    =  1+cosh(tau);
% 
n     =  2^12;  
t     = (0:2*pi/n:2*pi-2*pi/n).';
et    =    cosh(tau-i*t);
etp   = -i*sinh(tau-i*t);
%
zm =  0;
[zet,zetp,c] = mapdisk(et,etp,n,zm,'u');
%
Phio = fcau(et,etp,zet,zo,n,0);
%
nj    = 2^2;
nh    = n/nj;
th    = (2*pi/nh:2*pi/nh:pi-2*pi/nh).';
%% 
rv = [];  hv = [];
for kk=1:length(th)
    r = abs(cosh(tau-i*th(kk))-zo);
    %
    thet(1) =  -angle(zet(kk*nj+1));
    thet(2) =  -thet(1);
    %
    bet     =  -pi;
    sm = 0;
    for jj=1:2
        sm = sm+(-1)^(jj)*angle(i*(1+exp(-i*bet)*Phio)/(1-exp(-i*bet)*Phio)+cot((thet(jj)-bet)/2));
    end
    h  =  sm/pi;
    %
    rv = [rv;r];
    hv = [hv;h];
end
%
rv = [0  ;1 ;rv;2*cosh(tau)+1;5];
hv = [0  ;0 ;hv;1            ;1];
%%
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(rv,hv,'b','LineWidth',2)
hold on; box on
axis square
axis([0 5 0 1])
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[0:1:5]);
set(gca,'YTick',[0:0.2:1]);
xlabel('$r$','Interpreter','latex');
ylabel('$h(r)$','Interpreter','latex');
set(fig1,'PaperSize',[6  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_hfexellipse.pdf', '-dpdf', '-fillpage');
%%
