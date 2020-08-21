clc;clear;
load('D:\QTprojects\SV_QT\my_UDP_reception_f32\out_nlms_Frost_f32.txt')
%%
close all;

load('noise_voice34mic_4.mat','VarName1');
fs = 8000;
X = 80*VarName1/40000;
 
voic1 = X(1:length(X)/2);
voic2 = X(length(X)/2+1:end);
voic1(1:4000) = [];
voic2(1:4000) = [];

FF1 = fft(voic1);
FF2 = fft(voic2);
FF1(1) = 0;
FF2(1) = 0;
voic1 = ifft(FF1);
voic2 = ifft(FF2);

t_duration = length(voic1)/fs;
t = 0:1/fs:t_duration-1/fs;
f = (-length(voic1)/2:length(voic1)/2-1)*fs/length(voic1);




IN_BUFF_SIZE = 64;
OUT_BUFF_SIZE = 8;
T_INTERV = 20;
NUM_INTERV = 1000/T_INTERV;

s_in(:,1) = voic1(1:OUT_BUFF_SIZE*T_INTERV*NUM_INTERV);
s_in(:,2) = voic2(1:OUT_BUFF_SIZE*T_INTERV*NUM_INTERV);

Nel = 2;
c = 343;
ff = 2000;
d = c/ff/2;
J = 32;

out_MK = out_nlms_Frost_f32;

%% ...................................NLMS MK.................. ............
[y_nlms_MK, W_nlms_frost] = my_nlms_Frost_2sposob_2mic(s_in, J, 1);


%% 
figure()
hold on
plot(out_MK)
plot(y_nlms_MK)
% plot(y_nlms-y_nlms_MK)


%%
sound(out_MK,8e3)
%%
sound(y_nlms_MK,8e3)





