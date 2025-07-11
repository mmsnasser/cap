clear
addpath ../bie;addpath ../fmm;addpath ../files;
%%
n     =  2^10;
t     = (0:2*pi/n:2*pi-2*pi/n).';
%%
%
m     =  1;
alpha =  0.75i;
z1    =  0;
et    = []; etp =[];
et    =  exp(i*t);
etp   =  i*exp(i*t);
et    = [et  ;  0.75*cos(t)-0.5i*sin(t)];
etp   = [etp ; -0.75*sin(t)-0.5i*cos(t)];
% 
A     =  et-alpha;
%
gam   = -log(abs((et-z1)./(alpha-z1)));
[rho,h] = fbie(et,etp,A,gam,n,5,[],5e-14,50);
% 
for k=1:m+1
    Jk    = (k-1)*n+1:k*n;
    hk(k) = mean(h(Jk));
end
q  = exp(hk(2)-hk(1));
caph = q
%
%%