function y = mod2n(x,n)
if x>=0
 y = rem(x+2^(n-1),2^n)-2^(n-1);
else
 y = 2^(n-1)-1-rem(2^(n-1)-1-x,2^n);
end