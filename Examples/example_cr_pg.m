clear
addpath ../bie; addpath ../fmm; addpath ../files;
%%
% In this example, the function 
%
elv  =  [3,4,5,6,8,10,12,15,16,20,24,30,32].';
a    =   0.5;
% 
n        =   3*5*2^10;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
% 
for jj=1:length(elv)
%
el    =  elv(jj);
tht   =  2*pi/el;
for k=1:el
    ver(k,1)  = a*exp(-i*(k-1)*tht);
end
%
alphav   =  [ 0 ];
delt     =  [ 1 ];
alpha    =  0.5*(1+a)*exp(i*tht/2);
m        =   length(delt); 
%
% parametrization of \Gamma_1,\Gamma_2, the internal boundaries
eto   =    exp(i*t);
etop  =  i*exp(i*t);
% 
[eti,etip]=polygonp(ver,n/el);
%
et   = [eto ;eti ];
etp  = [etop;etip];
% 
% 
[capn(jj,1),~] = capg(et,etp,alphav,delt,alpha);
ra  =  a*sqrt(sin(tht))*sqrt(el/(2*pi));
rp  = (el*a/2)*sin(tht/2);
capa(jj,1) = 2*pi/log(1/ra);
capp(jj,1) = 2*pi/log(1/rp);
% 
%%
fig1=figure(1);clf;
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
if el==8
    print(fig1, 'fig_cr_pg.pdf', '-dpdf', '-fillpage');
end
% 
drawnow
%%
end

%%
fig2 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
plot(elv,capn,'-b','LineWidth',1.5)
hold on; box on
plot(elv,capa,':k','LineWidth',1.5)
plot(elv,capp,'-.r','LineWidth',1.5)
xlabel('$m$','FontSize',14,'Interpreter','latex');
% ylabel('Capacity','FontSize',14,'Interpreter','latex');
legend({'${\rm cap}(\Omega,E)$','${\rm cap}(\Omega,\overline{B^2(0,r_1)})$','${\rm cap}(\Omega,\overline{B^2(0,r_2)})$'},...
    'Location','east','Interpreter','latex')
% axis([1e1 1e4 1e-10  1e-0])
grid on; 
axis square
set(fig2,'PaperSize',[5  5]);
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(ax,'LooseInset',get(ax,'TightInset'))
%
print(fig2, 'fig_cirplg.pdf', '-dpdf', '-fillpage');%fillpage bestfit
% 
%%
erra = (capn-capa);
errp = (capn-capp);
fig3 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(elv,erra,'-b','LineWidth',1.5)
plot(elv,errp,'-r','LineWidth',1.5)
xlabel('$s$','FontSize',14,'Interpreter','latex');
% ylabel('Capacity','FontSize',14,'Interpreter','latex');
% axis([1e1 1e4 1e-10  1e-0])
grid on; 
axis square
set(fig3,'PaperSize',[5  5]);
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(ax,'LooseInset',get(ax,'TightInset'))
%
print(fig3, 'fig_cirplg2.pdf', '-dpdf', '-fillpage');
% 
%%