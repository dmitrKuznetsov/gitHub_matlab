%Обрати внимание на коэффициенты
clc;clear;close all;

load('data\micro12_noisAngle55.mat')

mic1(1:1000) = [];
mic2(1:1000) = [];
nFrost(1:1000) = [];




MAX = max([abs(mic1) abs(mic2) abs(nFrost)]);

% figure()
% plot(mic1/MAX)
% figure()
% plot(mic2/MAX)
% figure()
% plot(nFrost/MAX)
% 
figure(689)
hold on

plot(mic2/MAX)
plot(mic1/MAX)
plot(nFrost/MAX)
%%
% s_in(:,1) = mic1.';
% s_in(:,2) = mic2.';
% [y_nlms_MK, ~] = my_nlms_Frost_2sposob_2mic(s_in, 32, 0.1);
% figure(765)
% hold on
% 
% 
% plot(y_nlms_MK)
% plot(nFrost)


%%

% sound(mic1/MAX,8e3)
% %%
% sound(mic2/MAX,8e3)
% %%
% sound(nFrost/MAX,8e3)


sound(mic1/max(mic1),8e3)
%%
sound(mic2/max(mic2),8e3)
%%
sound(nFrost/max(nFrost),8e3)











