function   w = mobius(z,zv,wv)
% 
p  = @(x)((wv(3)-wv(1))*(wv(2)-wv(3))*(zv(2)-zv(1))*(x-zv(3)));
q  = @(x)((wv(2)-wv(1))*(zv(2)-zv(3))*(x-zv(1)));
r  = @(x)((wv(2)-wv(3))*(zv(2)-zv(1))*(x-zv(3)));
f  = @(x)(wv(3)+p(x)./(q(x)-r(x)));
%
w  =  f(z);
%
end