clc;clear;close all;
inPwr = 0.01; 
in = sqrt(inPwr)*randn(1000,1); 
mean(in.^2)

% System
order = 64;
LP = fir1(order-1,0.5,'low');
out = filter(LP, 1, in);


%% matlab RLS
lamda = 1;
rls = dsp.RLSFilter(order, 'ForgettingFactor' ,lamda);
[y_rls,e_rls] = rls(in,out);

%% my RLS
R_inv = 1000*diag(ones(order,1));
d = out;
alpha = zeros(length(d),1);
y = zeros(length(d),1);
x = zeros(order, 1);
h = zeros(order, 1);
for k = 1:length(d)
    x = [in(k); x(1:end-1)];
    g = R_inv*x/(   lamda + x'*R_inv*x   );
    R_inv = 1/lamda*(   R_inv - g*x'*R_inv   );
    y(k) = h'*x;
    alpha(k) = d(k) - y(k);
    h = h + g*alpha(k)';
end

figure(1)
hold on
plot(e_rls)
plot(alpha)

figure(2)
hold on
plot(y_rls)
plot(y)























