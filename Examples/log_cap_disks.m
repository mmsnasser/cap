clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n     =  2^10;
t     = (0:2*pi/n:2*pi-2*pi/n).';
%%
rv    = [linspace(2,8,61)].';
%%
for kk=1:length(rv)
r     =  rv(kk);
cent  =  [0    ; r    ;  r*i    ; -r   ; -r*i];
rad   =  [0.8  ; 1    ;  1      ;  1   ;  1  ];
%
m  =  length(rad)-1;
et = [];  etp = [];
for k=1:m+1 
    et  = [et  ; cent(k)+rad(k)*exp(-i*t)];
    etp = [etp ;      -i*rad(k)*exp(-i*t)];
end
% 
A     =  ones(size(et));
%
for k=1:m+1
    gam{k} = -log(abs((et)-cent(k)));
    [rho{k},h{k}] = fbie(et,etp,A,gam{k},n,5,[],5e-14,50);
end
% 
for k=1:m+1
    for j=1:m+1
        Jj = (j-1)*n+1:j*n;
        Mat(j,k) = mean(h{k}(Jj));
    end
end
Mat(1:m+1,m+2) = -1;
Mat(m+2,1:m+1) =  1;
Mat(m+2,m+2)   =  0;
delt  = zeros(m+2,1); delt(m+2) = 1;
sol   = Mat\delt; lkap = sol(m+2);
%
kap = exp(lkap)
%
logcap(kk,1) = kap;
end
%%

%%
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(rv,logcap,'b','LineWidth',2);
plot(rv,1.33+0.706*rv,':r','LineWidth',2);
axis square
% axis([-1  1 -1  1])
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick',[-3:1:3]);
% set(gca,'YTick',[-3:1:3]);
xlabel('$r$','Interpreter','latex');
% ylabel('${\rm cap}_l(E)$','Interpreter','latex');
legend({'${\rm cap}_l(E)$','$1.33+0.706r$'},'Interpreter','LaTeX',...
        'location','northwest');
set(fig1,'PaperSize',[5  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_logcap.pdf', '-dpdf', '-fillpage');
%%
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
plot(rv,1.33+0.706*rv-logcap,'k','LineWidth',2);
axis square
% axis([-1  1 -1  1])
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick',[-3:1:3]);
% set(gca,'YTick',[-3:1:3]);
xlabel('$r$','Interpreter','latex');
ylabel('$1.33+0.706r-{\rm cap}_l(E)$','Interpreter','latex');
set(fig1,'PaperSize',[5  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_logcapd.pdf', '-dpdf', '-fillpage');
%%