clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n     =  9*5*7*2^3;
t     = (0:2*pi/n:2*pi-2*pi/n).';
%%
zo    = 0; kv = [];  logcap =[];
for kk=3:10
    vert = sqrt(pi/(kk*sin(pi/kk)))*exp(-i*[0:1:kk-1]*2*pi/kk);
    [et,etp]=polygonp(vert,n/kk);
    %
    A     =  ones(size(et));
    gam   = -log(abs(et-zo));
    [~,h] =  fbie(et,etp,A,gam,n,5,[],5e-14,50);
    %
    Mat   = [mean(h) -1;1  0];
    delt  = [0;1];
    sol   =  Mat\delt; 
    lkap  =  sol(2);
    kap   =  exp(lkap);
    %
    kv      =  [kv;kk];
    logcap  =  [logcap;kap];
end
%%

%%
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(kv,logcap,'b','LineWidth',2);
% plot(kv,1+kv-kv,':r','LineWidth',2);
axis square
% axis([-1  1 -1  1])
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick',[-3:1:3]);
% set(gca,'YTick',[-3:1:3]);
xlabel('$k$','Interpreter','latex');
ylabel('${\rm cap}_l(E_k)$','Interpreter','latex');
% legend({'${\rm cap}_l(E)$','$1.2+0.74r$'},'Interpreter','LaTeX',...
%         'location','northwest');
set(fig1,'PaperSize',[5  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_logcap_Ep.pdf', '-dpdf', '-fillpage');
%%