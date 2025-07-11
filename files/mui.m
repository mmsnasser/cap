function x=mui(y)
% mui.m
% Compute the inverse normalized quotient of elliptic integrals;  mu^-1(r);
% where m(r)=(pi/2)*K'(r)/K(r)
% See Equation (5.1) in: G. D. Anderson, M. K. Vamanamurthy, and M.
% Vuorinen: Conformal invariants, inequalities and quasiconformal maps.  
% J. Wiley, 1997.
% See also page 126 in:
% P. Hariri, R. Klen, and M, Vuorinen, Conformal invariant metrics and
% quiasiconformal mappings, Springer, 2020.
%
% 
if y<.5
    x = sqrt(1-mui(pi^2/(4*y))^2);
else
    q    =  exp(-2*y);
    x    = (theta23(2,0,q,1e-16)/theta23(3,0,q,1e-16)).^2;
end
%
% p=1;
% for n=1:100
%     p=p*tanh((2*n-1)*y);
% end
% x = sqrt(1-p^8);
% 
% x=1;
% for n=1:100
%     x=x*tanh((2*n-1)*pi^2/(4*y))^4;
% end
%
end