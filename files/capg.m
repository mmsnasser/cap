function [cap,a,uz] = capg(et,etp,alphav,delt,alpha,z)
% Compute the capacity of the generalized condensers (G,F,delta) where G is
% a bounded domain
% Input:
% 1,2) et, etp: parametrization of the boundary and its first derivative
% 3) alphav=[alphav(1),...,alphav(mp)]: alphav(j) is an auxiliary point
% interior to \Gamma_j
% 4) deltv: the value of u on the boundary of E_k, k=1,2,...,m
% 5) alpha: alpha is an auxiliary point in G
% 6) zv: a vector of points in Omega = G\UE_k
% Output:
% cap , a , u(z)
%
% Computing the constants \h_{j,k} for j=1,2,...,m
%
m = length(delt);
n = length(et)/(m+1); 
A = et-alpha;
for k=1:m
    gamk{k}=log(abs(et-alphav(k)));
    [mu{k},h{k}]=fbie(et,etp,A,gamk{k},n,5,[],1e-14,100);
    for j=1:m+1
        hjk(j,k)=mean(h{k}(1+(j-1)*n:j*n));
    end
end
% Computing the constants a_k  for k=1,2,...,m
mat  =  hjk; mat(1:m+1,m+1)=1; 
rhs  = [0;delt(:)];
x    =  mat\rhs; a=x(1:m,1); c = x(m+1);
% Computing the capacity
cap  = (2*pi)*sum(delt.*a);
%
uz   = [];
%
if( nargin == 6 ) 
    Afet = zeros(size(et));
    for k=1:m
        Afet = Afet + a(k)*(gamk{k}+h{k}+i*mu{k});
    end
    fet = Afet./A;
    fzv = fcau(et,etp,fet,z);
    uz  = c+real((z-alpha).*fzv);
    for k=1:m
        uz = uz-a(k)*log(abs(z-alphav(k)));
    end
end
%
end