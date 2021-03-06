
%clc;clear;
%%
close all;

load('mic1.txt');


mic1 = reshape(mic1.',[],1);
mic1 = mic1(1:160);
order = 32;
% HP = fir1(order-1,0.5,'high');
load('ImResp_CompFIR.mat')

mic2 = (filter(bhi, 1, mic1));

nlms = dsp.LMSFilter('Length', order,'StepSize', 1,'Method','Normalized LMS','InitialConditions', zeros(32,1));

[out_matlab,e_matlab,w_matlab] = nlms(mic1,mic2);


Value_flip = flipud(Value1);
figure(36565)
hold on
% plot(LP/max(LP))
plot((w_matlab)/max(w_matlab))
plot(Value_flip/max(Value_flip))


% figure(6)
% freqz(Value_flip)




