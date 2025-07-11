clear
addpath ../bie; addpath ../fmm; addpath ../files;
%%
% In this example, the function 
%
n        =   2^11;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
r        =   0.1;
delt     =   ones(2,1);
%
alphavz  =  [0.5;-0.5];
alphaz   =   0.00;
% 
zo   =  alphavz(1);
T    =  @(z)((z-zo)./(1-conj(zo).*z));
Tp   =  @(z)((1-abs(zo)^2)./((1-conj(zo).*z).^2));
%
zet1   =  alphavz(1)+r*exp(-i*t);
zet1p  = -r*i*exp(-i*t);
% 
zet2   =  alphavz(2)+r*exp(-i*t);
zet2p  = -r*i*exp(-i*t);
% 
% parametrization of \Gamma_1,\Gamma_2, the internal boundaries
m     =   length(delt); 
eto   =    exp(i*t);
etop  =  i*exp(i*t);
%
et1   =   T(zet1);
et1p  =   Tp(zet1).*zet1p;
% 
et2   =   T(zet2);
et2p  =   Tp(zet2).*zet2p;
% 
et   = [eto  ;et1  ;et2  ];
etp  = [etop ;et1p ;et2p ];
% 
alpha  = T(alphaz);
alphav = T(alphavz);
%%
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=2:m+1
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'LineWidth',2)
end
k=1; crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
plot(real(crv),imag(crv),'b','LineWidth',2)
plot(real(alpha),imag(alpha),'pr','MarkerSize',8,'MarkerFaceColor','r')
plot(real(alphav),imag(alphav),'dr','MarkerSize',8,'MarkerFaceColor','r')
axis([-1.01  1.01 -1.01  1.01])
set(gca,'XTick',[-1:0.5:1]);
set(gca,'YTick',[-1:0.5:1]);
set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(ax,'LooseInset',get(ax,'TightInset'))
% print(fig1, 'fig_7d.pdf', '-dpdf', '-fillpage');
%%
% 
[capn,a,~] = capg(et,etp,alphav,delt,alpha);
capn
a
% 