% 
clear;clc
addpath ../bie;addpath ../fmm;addpath ../files;
%
n   =  2^13;
%
k     =  0.1059457511865234
b     =  1/mypsi(k)
%%
v     = [1     0    b*i    1+b*i];
%
modex =  mu(k)/pi;
%
n = 2^12
%
[et,etp]=polygonp(v,n/4);
%
zo = (v(1)+v(3))/2;
[zet,zetp] = mapdisk(et,etp,n,zo,'u');
vw = [zet(1) ; zet(n/4+1) ; zet(2*n/4+1) ; zet(3*n/4+1)];
modn = moddisk(vw(1),vw(2),vw(3),vw(4));
%
rerr = abs(modex-modn)/modex
%%
