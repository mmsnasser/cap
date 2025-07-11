clear
addpath ../bie; addpath ../fmm; addpath ../files;
%%
% In this example, the function 
%
n        =   2^13;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
a        =   0.5;
%
alphav   =  [ 0 ];
delt     =  [ 1];
alpha    =  (1+a)/2;
m        =   length(delt); 
% 
svo = a*[0.05:0.05:0.9,0.91:0.01:0.99].';
% 
for jj=1:length(svo)
% parametrization of \Gamma_1,\Gamma_2, the internal boundaries
eto   =   exp(i*t);
etop  = i*exp(i*t);
% 
s   =  svo(jj);
[c1,r1]=my3Pts( a,-s*i,-a);
[c2,r2]=my3Pts(-a, s*i, a);
ver      =  [ a  ; -a  ];
cent     =  [ c1 ;  c2 ];
dir      =  [-1  ; -1  ];
[eti,etip]=plgcirarcp(ver,cent,dir,n/2);
%
% 
et   = [eto ;eti ];
etp  = [etop;etip];
% 
% 
[capn(jj,1),~] = capg(et,etp,alphav,delt,alpha);
% 
estn(jj,1) = 2*pi/log(2*(pi-atan(a/(r1-s)))/(a*pi));
% 
end
%%
sv = [0;svo;a];
cap0 = 2*pi/mu(2*a/(1+a^2));  est0 = 2*pi/log(2/a);
cap1 = 2*pi/log(1/a);         est1 = 2*pi/log(1/a);
cap = [cap0;capn;cap1];
est = [est0;estn;est1];
%%
s   =  a/2;
[c1,r1]=my3Pts( a,-s*i,-a);
[c2,r2]=my3Pts(-a, s*i, a);
ver      =  [ a  ; -a  ];
cent     =  [ c1 ;  c2 ];
dir      =  [-1  ; -1  ];
[eti,~]=plgcirarcp(ver,cent,dir,n/2);
%
% 
et   = [eto ;eti ];
%
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
set(gca,'XTick',[-1:0.5:1]);
set(gca,'YTick',[-1:0.5:1]);
set(fig1,'PaperSize',[5  5]);
grid on; axis square
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
set(ax,'LooseInset',get(ax,'TightInset'))
%
print(fig1, 'fig_crlen.pdf', '-dpdf', '-fillpage');
% 
%%
fig2 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(sv,cap,'-b','LineWidth',1.5)
hold on; box on
plot(sv,est,':k','LineWidth',1.5)
xlabel('$s$','FontSize',14,'Interpreter','latex');
% ylabel('Capacity','FontSize',14,'Interpreter','latex');
legend({'${\rm cap}(\Omega,E)$','Estimiation'},...
    'Location','northwest','Interpreter','latex')
axis([0 0.5 4  10])
set(gca,'XTick',[0:0.1:0.5]);
grid on; 
axis square
set(fig2,'PaperSize',[5  5]);
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(ax,'LooseInset',get(ax,'TightInset'))
%
print(fig2, 'fig_cirlen.pdf', '-dpdf', '-fillpage');%fillpage bestfit
% 
%%
err = abs(cap-est);
fig3 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(sv,err,'-b','LineWidth',1.5)
hold on; box on
xlabel('$s$','FontSize',14,'Interpreter','latex');
% ylabel('Capacity','FontSize',14,'Interpreter','latex');
axis([0 0.5 0  0.03])
set(gca,'XTick',[0:0.1:0.5]);
grid on; 
axis square
set(fig3,'PaperSize',[5  5]);
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(ax,'LooseInset',get(ax,'TightInset'))
%
print(fig3, 'fig_cirlen2.pdf', '-dpdf', '-fillpage');
% 
%%