% clc;clear;
close all;

load('mic1.txt');


mic1 = reshape(mic1.',[],1);
in = mic1/2^(31-6);     %(1:160)
order = 32;

load('ImResp_CompFIR.mat')

d = filter(bhi, 1, in);


[out, e, b] = NLMS_f32_MK(in,d,order);


Value_flip = flipud(Value);
figure(1);
hold on
plot(bhi/max(bhi))
plot(b/max(b))
plot(Value_flip/max(Value_flip))

figure(2);
plot(e)

figure(3);
plot(out)


