function [out, e, b] = NLMS_f32_MK(in,d,order)
mu = 1;
e = zeros(length(d),1);
out = zeros(length(d),1);
x = zeros(order, 1);

b = zeros(order, 1);
for k = 1:length(d)
    x = [in(k); x(1:end-1)];
    out(k) = b'*x;
    
    e(k) = d(k) - out(k);
    E = (x' * x ) + 0.000000119209289;
    b = b + mu / E * x * e(k)';
end


end

