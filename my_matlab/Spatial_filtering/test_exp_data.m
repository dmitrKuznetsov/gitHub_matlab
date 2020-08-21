%Обрати внимание на коэффициенты
clc;clear;close all;

load('micro12_noisAngle55.mat')

% mic1(1:1000) = [];
% mic2(1:1000) = [];
% nFrost(1:1000) = [];




MAX = max([abs(mic1) abs(mic2) abs(nFrost)]);



% figure(689)
% hold on
% plot(mic2/MAX)
% plot(mic1/MAX)
% plot(nFrost/MAX)
%%
s_in(:,1) = mic1.'/100;
s_in(:,2) = mic2.'/100;
J = 200;
[y_nfrost, W_nfrost] = func_LC_NLMS_2sposob_2mic_imitac_single(s_in, J, 0.1); %func_Frost(s_in, 32, 2, 0.005);  %func_LC_NLMS_2sposob_2mic_imitac_single(s_in, J, 0.1);
% [y_lcRLS, W_lcRLS] = func_LC_RLS(s_in, 2, J, 1);
W_SMI = func_SMI(mic1,mic2,J);
y_SMI = apply_coef(s_in, W_SMI, J, 2);
Nel = 2;
figure(765)
hold on
% plot(y_lcRLS)
plot(y_nfrost)
plot(y_SMI)

plotBP_from_freq(W_SMI, 2, J, 0, -55);

% plotSpectr(mic2, 1024, 7667, 8e3)
%%


% sound(mic1/max(mic1),8e3)
% %%
% sound(mic2/max(mic2),8e3)
% %%
% sound(nFrost/max(nFrost),8e3)
%%
sound(y_nfrost(1000:end)/max(y_nfrost(1000:end)),8e3)

%%
sound(y_SMI(1000:end)/max(y_SMI(1000:end)),8e3)







