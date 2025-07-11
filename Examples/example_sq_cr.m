clear
addpath ../bie; addpath ../fmm; addpath ../files;
%%
% In this example, the function 
%
nv = [8:4:100,110:10:200,220:20:400,440:40:1000].';
% 
for jj=1:length(nv)
n        =   nv(jj);
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
a        =   0.5;
%
alphav   =  [ 0 ];
delt     =  [ 1];
alpha    =  (1+a)/2;
% 
m        =   length(delt); 
%
% parametrization of \Gamma_1,\Gamma_2, the internal boundaries
vero     =  [ 1+i; -1+i ; -1-i ; 1-i];
cento    =  [ inf;  inf ;  inf ; inf];
diro     =  [ 0  ;  0   ;  0   ; 0  ];
[eto,etop]=plgsegcirarcp(vero,cento,diro,n/4);
%
eti   =    a*exp(-i*t);
etip  = -i*a*exp(-i*t);
% 
et   = [eto ;eti ];
etp  = [etop;etip];
% 
% 
[capn,~] = capg(et,etp,alphav,delt,alpha);
% 
end
%%
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:m+1
    crv    =  et((k-1)*n+1:k*n,1); crv(n+1)  =  crv(1);
    plot(real(crv),imag(crv),'b','LineWidth',2)
end
plot(real(alpha),imag(alpha),'pr','MarkerFaceColor','r','MarkerSize',12)
plot(real(alphav),imag(alphav),'dk','MarkerFaceColor','k','MarkerSize',10)
axis([-1.01  1.01 -1.01  1.01])
set(gca,'XTick',[-1:0.5:1],'FontSize',14);
set(gca,'YTick',[-1:0.5:1]);
set(fig1,'PaperSize',[10 10]);
grid on; axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(ax,'LooseInset',get(ax,'TightInset'))
%
print(fig1, 'fig_sq_cr.pdf', '-dpdf', '-fillpage');
% 
%%
