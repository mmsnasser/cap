clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
rv    =  [linspace(0.2,5,30)]';
nv    =  [2^8,2^10,2^12];
%%
for jj=1:length(nv)
    n   =  nv(jj);
    t   = (0:2*pi/n:2*pi-2*pi/n).';
for kk=1:length(rv)
    r = rv(kk);
et    =    cosh(r-i*t);
etp   = -i*sinh(r-i*t);
%
%
zo =  0;
[zet,zetp,c] = mapdisk(et,etp,n,zo,'u');
redmodn =  -log(c)/(2*pi);
% 
redmode =  (log(2)-r)/(2*pi);
%
rerrkk = abs((redmode-redmodn)/redmode);
rerrkk(rerrkk<eps)=eps;
rerr(kk,jj) = rerrkk;
end
end
%%
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
semilogy(rv,rerr(:,1),'k','LineWidth',2)
hold on; box on
semilogy(rv,rerr(:,2),'b','LineWidth',2)
semilogy(rv,rerr(:,3),'r','LineWidth',2)
axis square
axis([0 5 1e-16 1e-12])
set(gca,'LooseInset',get(gca,'TightInset'))
set(gca,'XTick',[0:1:5]);
set(gca,'YTick',[1e-16,1e-15,1e-14,1e-13,1e-12]);
xlabel('$r$','Interpreter','latex');
ylabel('Relative Error','Interpreter','latex');
legend({'$n=2^8$','$n=2^{10}$','$n=2^{12}$'},'Interpreter','LaTeX',...
        'location','northeast');
set(fig1,'PaperSize',[6  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_rerrexellipse.pdf', '-dpdf', '-fillpage');
%%
