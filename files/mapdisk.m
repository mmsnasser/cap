function [zet,zetp,c] = mapdisk(et,etp,n,zz,type)
% mapdisk.m
% 9-5-2021
% 1) type=='b' if G is bounded and zz=alpha is a given point in G:
% the function mapdisk(et,etp,n,zz,type) computes the conformal mapping
% w = Phi(z) from the domain G onto the unit disk with the normalizations
% Phi(alpha)=0 and c=Phi'(alpha)>0.
% 2) type=='u' if G is unbounded and zz=z1 is exterior to G:
% the function mapdisk(et,etp,n,zz,type) computes the conformal mapping
% w = Phi(z) from the unbounded domain onto the unit disk with
% the normalizations Phi(inf)=0 and c=Phi'(inf)>0. 
%
% Input:
% 1) et, etp: the parametrization of the boundary of G and its derivative
% 2) zz: zz=alpha is a given point G if G is bounded and zz=zo is a point 
%           in the exterior of G if G is unbounded
% 
% Output:
% 1) zet=Phi(et) and zetp is the derivative of zet.
% 2) c=Phi'(alpha)>0 for bounded G and c=Phi'(inf)>0 for unbounded G
% 
% Note that the boundary of G must be oriented such that G is on the 
% left of its boundary.
% 
%
if type=='b' 
    A = et-zz; k = 1;
elseif type=='u' 
    A = ones(size(et)); k = -1;
end
gam     = -k*log(abs(et-zz));
% 
[rho,h] =  fbie(et,etp,A,gam,n,5,[],1e-14,50);
c       =  exp(-mean(h));
Af      =  gam+mean(h)+i*rho;
if type=='b' 
    zet = c*(et-zz).*exp(Af);
elseif type=='u' 
    zet = (c./(et-zz)).*exp(Af);
end
zetp   =   derfft(real(zet))+i*derfft(imag(zet));
%%
end