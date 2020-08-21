clc;clear;close all;


inPwr = 0.01; 
in = sqrt(inPwr)*randn(1000,1); 
% mean(in.^2)

% System
order = 64;
LP = fir1(order-1,0.5,'low');
out = filter(LP, 1, in);


%% matlab NLMS
mu = 0.5;
lms = dsp.LMSFilter(order,'StepSize',mu,'Method', 'Normalized LMS');
[y_lms,e_lms,w_lms] = lms(in, out);

% figure(1)
% hold on
% stem(LP)
% stem(w_lms)
% 
% figure(2)
% plot(e_lms)
%% my NLMS

beta = 2.2204460492503131e-016;
d = out;
alpha = zeros(length(d),1);
y = zeros(length(d),1);
x = zeros(order, 1);
h = zeros(order, 1);
for k = 1:length(d)
    x = [in(k); x(1:end-1)];
    y(k) = h'*x;
    alpha(k) = d(k) - y(k);
    h = h + mu / (x' * x + beta) * x * alpha(k)';
    
end

% figure(3)
% hold on
% stem(LP)
% stem(h)
% 
% figure(4)
% plot(alpha)


%% Сравнение
figure(5)
hold on
plot(alpha)
plot(e_lms)

figure(6)
hold on
plot(y_lms)
plot(y)


























