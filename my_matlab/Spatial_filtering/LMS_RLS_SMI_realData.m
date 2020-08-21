clc;clear;close all;
load('noise_voice34mic_4.mat','VarName1');
fs = 8000;
X = 80*VarName1/40000;
% X(1:300) = []; 
% voic1 = X; 
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

Nel = 2;
c = 343;
ff = 2000;
d = c/ff/2;
hap = dsp.AudioPlayer('SampleRate',fs);


s_in(:,1) = voic1;
s_in(:,2) = voic2;




%%
% len_posled = 50:5:300;
len_posled = 1000;
% LMS 
L = 32;
K = Nel;
mu = 0.05;

for lp = 1: length(len_posled)
    lp

    %LMS
    [y_LMS,W_LMS] = my_Frost(s_in(1:len_posled(lp),:), L, K, mu);
    %RLS 
    [y_LC_RLS, W_RLS] = func_LC_RLS_broadband(s_in(1:len_posled(lp),:),K,L,1);
    
    %SMI
    W_SMI = SMI_func(s_in(1:len_posled(lp),1),s_in(1:len_posled(lp),2),L);
    
    
    [out_LMS] = apply_coef(s_in, W_LMS, L, K);
    [out_RLS] = apply_coef(s_in, W_RLS, L, K);
    [out_SMI] = apply_coef(s_in, W_SMI, L, K);


end
%%
figure(12112)
hold on
plot(out_LMS)
plot(out_RLS)
plot(out_SMI)
legend('LMS','RLS','SMI')


audiowrite('raw_averaged.wav',mean(s_in,2),fs)
audiowrite('out_LMS.wav',out_LMS,fs)
audiowrite('out_RLS.wav',out_RLS,fs)
audiowrite('out_SMI.wav',out_SMI,fs)

%%
figure(1);
plot(t, voic1);
hold on
plot(t, voic2);
ylim([-0.8 0.8])