clc;clear;close all;
load('D:\QTprojects\SV_QT\my_UDP_reception_f32\out_nlms_Frost_f32.txt')
sig = out_nlms_Frost_f32(1000:end)/max(abs(out_nlms_Frost_f32(1000:end)));
plot(sig)







sound(sig,8e3)









