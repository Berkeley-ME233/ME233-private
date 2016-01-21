function [R,S] = bezout(A,B,C)
% function [R,S] = bezout(A,B,C)
%
% solution of the Bezout Equation
% C(q-1) = A(q^(-1)) * R(q^(-1)) + q^(-1) * B(q^(-1))* S(q^(-1))
%
% A, B and C are polynomials in  q^(-1) = 1/q
% A and C must be polynomials with 1 as their constant coefficient
%       (i.e. A(1) = 1 and C(1) = 1)
% A and B must be co-prime
%
% returns polynomials R and S
if (nargin < 1)
    error('No enough arguments')
end
a = unpad(A);
b = unpad(B) ; 
c = unpad(C);

if c(1) ~= 1 || a(1) ~= 1
    error('either a or c do not have 1 as the constant coefficient')
end
n = length(a)-1; m=length(b)-1; nc = length(c)-1;
nr = m;
ns=max([n-1,nc-m-1]);
nv = m + ns + 1;
ce = [c(2:end) zeros(1,nv - nc )];
ae = [a(2:end) zeros(1,nv - n )];
D=zeros(nv,nv);
for jj=1:m
    D(jj:jj+n,jj) = a';
end
j=1;
for jj=m+1:nv
    D(j:j+m,jj) = b';
    j=j+1;
end
v = (D \ (ce - ae)')';
R = [1 v(1:m)];
S = v(m+1:end);

    


function vout = unpad(v)
if (size(v,1) > 1)
    v=v';
end
if (size(v,1) > 1)
    error('not a vector')
end
r1 =roots(v);
num_zeros = sum(r1==0);
vout = v(1:end-num_zeros);

