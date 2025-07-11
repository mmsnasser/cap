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
eto   =    exp(i*t);
etop  =  i*exp(i*t);
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
cape = 2*pi/log(1/a);
% 
rerr(jj,1) = abs(capn-cape)/cape;
% 
fprintf('%20d %20.14f %20.14f %20.4e\n',[n capn cape rerr(jj)])
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
print(fig1, 'fig_2cr.pdf', '-dpdf', '-fillpage');
% 
%%
fig2 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
loglog(nv,rerr,'-b','LineWidth',1.5)
hold on; box on
loglog(nv,1.5*exp(-0.295*nv),'-.k','LineWidth',1.5)
xlabel('$n$','FontSize',14,'Interpreter','latex');
ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'Relative Error','$O(e^{-0.295n})$'},...
    'Location','northeast','Interpreter','latex')
axis([1e1 1e3 1e-17  1e-0])
grid on; 
axis square
set(fig2,'PaperSize',[5 5]);
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))
%
print(fig2, 'fig_cirerr.pdf', '-dpdf', '-fillpage');%fillpage bestfit
% 
%%
