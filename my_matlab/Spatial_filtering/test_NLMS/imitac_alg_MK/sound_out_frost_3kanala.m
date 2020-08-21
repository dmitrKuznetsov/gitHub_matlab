%������ �������� �� ������������
clc;clear;close all;
% load('D:\QTprojects\SV_QT\my_UDP_reception_f32\out_nlms_Frost_f32.txt')
load('data\out_nlms_Frost_f32.txt')
OUT_BUFF_SIZE = 8;
T_INTERV = 20;
NUM_INTERV = 15000/T_INTERV;

PCM_reshaped = reshape(out_nlms_Frost_f32.', T_INTERV*OUT_BUFF_SIZE, []);
resh_mic1 = PCM_reshaped(:,1:4:end);
resh_mic2 = PCM_reshaped(:,2:4:end);
resh_nFrost = PCM_reshaped(:,3:4:end);


mic1 = reshape(resh_mic1, 1, []);
mic2 = reshape(resh_mic2, 1, []);
nFrost = reshape(resh_nFrost, 1, []);

mic1(1:1000) = [];
mic2(1:1000) = [];
nFrost(1:1000) = [];

% mic1(32001:end) = [];
% mic2(32001:end) = [];
% nFrost(32001:end) = [];


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











