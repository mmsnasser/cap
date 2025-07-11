clear
addpath ../bie; addpath ../fmm; addpath ../files;
%%
% In this example, the function 
%
n        =   2^11;
t        =  (0:2*pi/n:2*pi-2*pi/n).';
%
rv = [0.01:0.0025:0.17]; av = [];
for kk=1:length(rv)
r        =   rv(kk);
%
for k=1:7
    alphav(k,1)   =  (0.1+k/10)*exp(i*(k-1)*pi/2);
end
delt     =   ones(7,1);
alpha    =  -0.05;
% 
m        =   length(delt); 
%
% parametrization of \Gamma_1,\Gamma_2, the internal boundaries
et   =    exp(i*t);
etp  =  i*exp(i*t);
%
for k=1:7
    eti   =  alphav(k)+r*exp(-i*t);
    etip  = -r*i*exp(-i*t);% 
    % 
    et   = [et ;eti ];
    etp  = [etp;etip];
end
% 
%%
if r==0.1
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
print(fig1, 'fig_7d.pdf', '-dpdf', '-fillpage');
end 
%%
% 
[capn(kk,1),a,~] = capg(et,etp,alphav,delt,alpha);
av = [av a];
% 
end
% 
%%
% 
%%
fig2 = figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');  
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
for k=1:7
    plot(rv,av(k,:),'LineWidth',1.5)
end
xlabel('$r$','FontSize',14,'Interpreter','latex');
% ylabel('Relative Error','FontSize',14,'Interpreter','latex');
legend({'$a_1$','$a_2$','$a_3$','$a_4$','$a_5$','$a_6$','$a_7$'},...
    'Location','northwest','Interpreter','latex')
axis([0.05 0.17 0  2])
grid on; 
axis square
set(fig2,'PaperSize',[5 5]);
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',16)
set(gca,'LooseInset',get(gca,'TightInset'))
print(fig2, 'fig_7dak.pdf', '-dpdf', '-fillpage');%fillpage bestfit 
%%
