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
load('mic1.txt');
load('mic2.txt');

mic1 = reshape(mic1.',[],1);
mic2 = reshape(mic2.',[],1);
nlms = dsp.LMSFilter('Length', 32,'StepSize', 1,'Method','Normalized LMS','InitialConditions', ones(32,1));

[out_matlab,e_matlab,~] = nlms(mic1,mic2);              %(x,d)






figure(1)
hold on
plot(out_matlab/max(out_matlab)*2)
plot(out_nlms/max(out_nlms))
ylim([-1 1])

figure(2)
hold on

plot(e_matlab/max(e_matlab)*15)
plot(err_nlms/max(err_nlms))
ylim([-1 1])



figure(546)
hold on
plot(mic2.'-(out_nlms+err_nlms))
% plot(out_nlms+err_nlms)





