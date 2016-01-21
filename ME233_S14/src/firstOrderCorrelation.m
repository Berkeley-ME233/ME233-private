t = -8:1:8;
a1 = 0.8;
a2 = 0.7;
a3 = 0.6;
A = zeros(length(t),4);
A(:,1) = t';
for ii = -8:8
    A(ii+9,2) = a1^abs(ii);
    A(ii+9,3) = a2^abs(ii);
    A(ii+9,4) = a3^abs(ii);
end
A
