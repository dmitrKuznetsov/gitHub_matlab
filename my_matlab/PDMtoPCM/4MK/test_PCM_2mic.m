clc;clear;close all;
load('D:\QTprojects\SV_QT\my_UDP_v1\PCM_2m.txt')
out_buf_size = 8;
t_interv = 20;
num_interv = 5000/t_interv;       %!!!

(1:8000)
for ini = 0: (num_interv-1)
    mic1(1+ini*out_buf_size*t_interv   :  (ini+1)*out_buf_size*t_interv) ...
  = PCM_2m(1+2*ini*out_buf_size*t_interv   :  (2*ini+1)*out_buf_size*t_interv);
    
    mic2(1+ini*out_buf_size*t_interv    :   (ini+1)*out_buf_size*t_interv) ...
  = PCM_2m(1+(2*ini+1)*out_buf_size*t_interv   :  (2*ini+2)*out_buf_size*t_interv);
end

figure(1)
hold on

% ylim([-1 1])
% figure(2)
plot(mic1/max(mic1))
plot(mic2/max(mic1))
ylim([-1 1])
%%
sound(mic1/max(mic1),8000) %
%%
sound(mic2/max(mic2),8000) %
