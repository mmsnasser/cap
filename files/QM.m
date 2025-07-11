function md = QM(A,B)
% 
%
%
beta = @(x,y)(gamma(x)*gamma(y)/gamma(x+y));
%      
a =  1-(angle(A-1)-angle(A-B))/pi; b =  angle(B)/pi;
c = (pi-angle(A-1)+angle(B))/pi;
L = (beta(c-b,1-a)/beta(b,c-b))*exp(i*(b+1-c)*pi);
%
f = @(x)(h(x)); ff= @(x)(h(x)^2);
if f(1e-6)*f(1-1e-13)<0
    r = fzero(f,[1e-6,1-1e-13]);
else
    r = fminbnd(ff,0,1,optimset('TolX',1e-14));
end
md = (2/pi)*mu(r);
%
% 
function y = h(x)
   nn = ((1-x^2)^(c-a-b))*hypergeom([c-a,c-b],c+1-a-b,1-x^2);
   dd =  hypergeom([a,b],c,x^2);
    y =  nn/dd-real((A-1)/L);
end
end