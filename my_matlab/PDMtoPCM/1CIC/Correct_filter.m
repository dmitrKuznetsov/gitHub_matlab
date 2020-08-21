clc;clear;close all;
R = 32;          %коэф децим
D = 2;          %задержка
N = 4;          %порядок фильтра

fs = 8000*2;



in = [0 1 ones(1,16*1200-2)];
out_CIC_nmy = CIC_nmy_N4(in,R,37);

% figure(555)
% hold on
% plot(out_CIC_nmy)
% plot(out_CIC_mat)

H_my = abs(fft(out_CIC_nmy));


H_my(1) = 0;
H_my = H_my/max(H_my);
H_my(1) = 1;


H_ideal = [ones(1,length(H_my)/4 ) zeros(1,length(H_my)/2 ) ones(1, length(H_my)/4)];
% figure(666)
% hold on
% plot(H_my)
% plot(H_matlab)
% plot(H_ideal)



H_correct = H_ideal./H_my;
H_correct(length(H_correct)/2+1:end) = [];
figure(777)
plot((H_correct))
% 
% h_correct = ifft(H_correct,256);
% figure(888)
% plot(H_correct)


f = (0:length(H_correct))/length(H_correct);

bhi = fir2(256-1,f,[H_correct 0] );                                      %!!!!!
figure(888)
freqz(bhi,1)
F_resp_cor = abs(fft(bhi));
F_resp_cor_dB = 20*log10(F_resp_cor/max(F_resp_cor));



F_plot_real = F_resp_cor_dB(1:length(F_resp_cor_dB)/2);
f_sootv = (1:length(F_plot_real))/length(F_plot_real)*length(H_correct);


figure(999)
hold on
plot(f_sootv,F_plot_real)
plot(20*log10(H_correct/max(H_correct)))













