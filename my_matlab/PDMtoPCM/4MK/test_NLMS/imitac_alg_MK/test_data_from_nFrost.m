% данные хранятся на мк 
clc;clear;close all;
load('out_nlms_Frost_f32_data.mat','out_nlms_Frost_f32')
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

mic1(1:24000) = [];
mic2(1:24000) = [];
nFrost(1:24000) = [];
mic1(8001:end) = [];
mic2(8001:end) = [];
nFrost(8001:end) = [];

MAX = max([abs(mic1) abs(mic2) abs(nFrost)]);

s_in(:,1) = mic1.';
s_in(:,2) = mic2.';
[y_nlms_MK_imitac, ~] = my_nlms_Frost_2sposob_2mic_imitac_single(s_in, 32, 1);

load('D:\QTprojects\SV_QT\my_UDP_reception_f32\out_nlms_Frost_f32.txt');


figure(765)
hold on

plot(out_nlms_Frost_f32)
plot(y_nlms_MK_imitac)





