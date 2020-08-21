clc;clear;close all;
load('D:\QTprojects\SV_QT\my_UDP_v1\PCM_2m.txt')
out_buf_size = 8;
t_interv = 20;
num_interv = 1000/t_interv;       %!!!


for ini = 0: (num_interv-1)
    out_nlms(1+ini*out_buf_size*t_interv   :  (ini+1)*out_buf_size*t_interv) ...
  = PCM_2m(1+2*ini*out_buf_size*t_interv   :  (2*ini+1)*out_buf_size*t_interv);
    
    err_nlms(1+ini*out_buf_size*t_interv    :   (ini+1)*out_buf_size*t_interv) ...
  = PCM_2m(1+(2*ini+1)*out_buf_size*t_interv   :  (2*ini+2)*out_buf_size*t_interv);
end



%%
% sound(out_nlms/max(out_nlms),8000) %
%% matlab
close all;
load('mic1.txt');
% load('mic2.txt');

mic1 = reshape(mic1.',[],1);
% mic2 = reshape(mic2.',[],1);
order = 32;
LP = fir1(order-1,0.5,'high');

mic2 = floor(filter(LP, 1, mic1)*10^4);

nlms = dsp.LMSFilter('Length', order,'StepSize', 1,'Method','Normalized LMS','InitialConditions', zeros(32,1));

[out_matlab,e_matlab,w_matlab] = nlms(mic1,mic2);






figure(1)
hold on

 plot(out_nlms/max(out_nlms))
 plot(out_matlab/max(out_matlab))
% ylim([-1 1])

figure(2)
hold on

plot(err_nlms/max(err_nlms))
plot(e_matlab/max(e_matlab))

% ylim([-1 1])

% 
% figure(3)
% % plot(fftshift(w_matlab))
% % plot(LP)
% freqz(w_matlab)
% 
% figure(4)
% freqz(LP)
% 
% figure(5)
% freqz(LP)
% 
% Value_lr = flipud(Value);
% figure(676767)
% hold on
% plot(Value_lr)
% plot(Value)
% 
% figure(6)
% freqz(Value_lr)
% 
% figure(36565)
% hold on
% plot((w_matlab)/max(w_matlab))
% plot(LP/max(LP))
% plot(Value_lr/max(Value_lr))
% % 
% 
