function [zet,zetp,c] = diskcirslitmap(et,etp,n,alpha)
% ???.m
% 13-6-2021
% 
%
m       =  length(et)/n-1;
A       =  et-alpha; 
gam     = -log(abs(et-alpha));
% 
[rho,h] =  fbie(et,etp,A,gam,n,5,[],1e-14,50);
for k=1:m+1
    Jk = (k-1)*n+1:k*n;
    hk(k) = mean(h(Jk));
end
c       =  exp(-hk(1));
Af      =  gam+mean(h)+i*rho;
zet = c*(et-alpha).*exp(Af);
zetp   =   derfft(real(zet))+i*derfft(imag(zet));
%%
end