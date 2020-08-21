clc;clear;close all;
load('D:\QTprojects\SV_QT\my_UDP_v1\PCM_2m.txt')
out_buf_size = 8;
t_interv = 20;
num_interv = 1000/t_interv;       %!!!


for ini = 0: (num_interv-1)
    out_fir1(1+ini*out_buf_size*t_interv   :  (ini+1)*out_buf_size*t_interv) ...
  = PCM_2m(1+2*ini*out_buf_size*t_interv   :  (2*ini+1)*out_buf_size*t_interv);
    
    out_fir2(1+ini*out_buf_size*t_interv    :   (ini+1)*out_buf_size*t_interv) ...
  = PCM_2m(1+(2*ini+1)*out_buf_size*t_interv   :  (2*ini+2)*out_buf_size*t_interv);
end



load('mic1.txt');


mic1 = reshape(mic1.',[],1);

order = 32;
% HP = fir1(order-1,0.5,'high');
load('ImResp_CompFIR.mat')

mic2 = (filter(bhi, 1, mic1));

figure()
hold on

plot(mic2/max(mic2))
plot(out_fir1/max(out_fir1))