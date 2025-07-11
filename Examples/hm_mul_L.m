clc; clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n     =  2^11;
t     = (0:2*pi/n:2*pi-2*pi/n).';
%
cent  =  [0     ; 0.6  ; -0.3+0.5i ; -0.3-0.5i];
rad   =  [0.15  ; 0.2  ;  0.2      ;  0.2     ];
%%
m  =  length(rad);
et = [];  etp = [];
et    =    exp(i*t);
etp   =  i*exp(i*t);
for k=1:m 
    et  = [et  ; cent(k)+rad(k)*exp(-i*t)];
    etp = [etp ;      -i*rad(k)*exp(-i*t)];
end
% 
alpha =  0.5+0.5i;
A     =  et-alpha;
%%
for k=1:m
    gam{k} = log(abs((et)-cent(k)));
    [rho{k},h{k}] = fbie(et,etp,A,gam{k},n,5,[],5e-14,50);
end
% 
for k=1:m
    for j=1:m+1
        Jj = (j-1)*n+1:j*n;
        Mat(j,k) = mean(h{k}(Jj));
    end
end
Mat(1:m+1,m+1) = 1;
delt1 = zeros(m+1,1); delt1(2) = 1;
delt2 = zeros(m+1,1); delt2(3) = 1;
sol1 = Mat\delt1; a1 = sol1(1:m); c1 = sol1(m+1);
sol2 = Mat\delt2; a2 = sol2(1:m); c2 = sol2(m+1);
g1et = zeros(size(et));
g2et = zeros(size(et));
for k=1:m
    g1et = g1et+a1(k)*(gam{k}+h{k}+i*rho{k})./A;
    g2et = g2et+a2(k)*(gam{k}+h{k}+i*rho{k})./A;
end
%%
x = linspace(-1,1,251); 
y = linspace(-1,1,251); 
[xx,yy] = meshgrid(x,y); 
zz = xx+1i*yy;
zv = zz(:);
zv(abs(zv)>=1) = NaN+i*NaN;
for k=1:m
    zv(abs(zv-cent(k))<=rad(k)) = NaN+i*NaN;
end
z  = zv(abs(zv)>=0); z=z(:).';
%%
g1z = fcau(et,etp,g1et,z);
g2z = fcau(et,etp,g2et,z);
%
sig1z = real((z-alpha).*g1z)+c1;
sig2z = real((z-alpha).*g2z)+c2;
for k=1:m
    sig1z = sig1z-a1(k)*log(abs(z-cent(k)));
    sig2z = sig2z-a2(k)*log(abs(z-cent(k)));
end
% 
sig1v = NaN(size(zv));
sig1v(abs(zv)>=0)=sig1z;
sig1  = reshape(sig1v,size(zz));
% 
sig2v = NaN(size(zv));
sig2v(abs(zv)>=0)=sig2z;
sig2  = reshape(sig2v,size(zz));
%%
cntv  =  [0,0.01,0.05:0.1:0.95,1];
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
contourf(xx,yy,sig1,cntv,'color','k','linewidth',1);
colormap jet
colr = 0.8.*colormap+0.2.*ones(size(colormap));
brighten(colr,0.4)
caxis([0  1])
for k=1:m+1
    Jk = (k-1)*n+1:k*n;
    plot(real(et(Jk)),imag(et(Jk)),'b','LineWidth',2);
end
axis equal
axis([-1  1 -1  1])
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick',[-3:1:3]);
% set(gca,'YTick',[-3:1:3]);
set(fig1,'PaperSize',[5  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_hmL1.pdf', '-dpdf', '-fillpage');
%%
cntv  =  [0,0.01,0.05:0.1:0.95,1];
fig1=figure;
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
hold on; box on
contourf(xx,yy,sig2,cntv,'color','k','linewidth',1);
colormap jet
colr = 0.8.*colormap+0.2.*ones(size(colormap));
brighten(colr,0.4)
caxis([0  1])
for k=1:m+1
    Jk = (k-1)*n+1:k*n;
    plot(real(et(Jk)),imag(et(Jk)),'b','LineWidth',2);
end
axis equal
axis([-1  1 -1  1])
set(gca,'LooseInset',get(gca,'TightInset'))
% set(gca,'XTick',[-3:1:3]);
% set(gca,'YTick',[-3:1:3]);
set(fig1,'PaperSize',[5  5]);
grid on; grid('minor')
ax=gca; 
set(ax,'xminorgrid','on','yminorgrid','on')
ax.GridAlpha=0.25; ax.MinorGridAlpha=0.25;
set(gca,'FontSize',20)
%
print(fig1, 'fig_hmL2.pdf', '-dpdf', '-fillpage');
%%