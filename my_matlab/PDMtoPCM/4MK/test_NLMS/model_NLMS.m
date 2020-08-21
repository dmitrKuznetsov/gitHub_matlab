% clc;clear;
close all;

load('mic1.txt');
mic1 = reshape(mic1.',[],1);
in = mic1;      %(1:32)

load('ImResp_CompFIR.mat')
d = floor(filter(bhi, 1, mic1)/(2^31));

order = 32;
Value_flip = flipud(Value);

mu = 2^31-1;
%  mu = 1;
e = zeros(length(d),1);
y = zeros(length(d),1);
x = zeros(order, 1);
% b = zeros(order, 1);
b = ones(order, 1);
for k = 1:length(d)
    x = [in(k); x(1:end-1)];
    y(k) = floor(b'*x/(2^31));
    
    e(k) = d(k) - y(k);
%     E = mod2n((x' * x ),32);
    E = (x' * x )/(2^31);
    b = b + floor(mu / E * x * e(k)'/(2^31));
end


figure();
hold on
%plot(bhi/max(bhi))
plot(b/max(b))
plot(Value/max(Value))


%%
%     y_filt = (filter(b, 1, x))   
%     b'*x
%     y(k) = y_filt(k)

%     b.'



































