clc;clear;close all;


inPwr = 0.01; 
in = sqrt(inPwr)*randn(1000,1); 
mean(in.^2)

% System
order = 64;
LP = fir1(order-1,0.5,'low');
out = filter(LP, 1, in);


% LMS
mu = 1;
lms = dsp.LMSFilter(order,'StepSize',mu);
[y_lms,e_lms,w_lms] = lms(in, out);

figure(1)
hold on
stem(LP)
stem(w_lms)

figure(2)
plot(e_lms)

%RLS
rls = dsp.RLSFilter(order);
[y_rls,e_rls] = rls(in,out);

figure(4)
plot(e_rls)









