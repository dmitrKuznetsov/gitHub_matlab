% X = VarName1;
% PCM = X(1:16000);
% PDM = X(16001:end);
% save('PCM_raw.mat','PCM')
% save('PDM_raw.mat','PDM')
%%
clc; clear; close all;
load('PCM_raw.mat','PCM')
load('PDM_raw.mat','PDM')
pcm_err_ms = 8*50; % первые 50 мс отрезаем 
pdm_err_ms = 64*50;
PCM(1:pcm_err_ms) = [];
PDM(1:pdm_err_ms) = [];

% PDM_bi =  de2bi(PDM) ;
PDM_bi = fliplr( de2bi(PDM) );

PDM_bi_reshape = reshape(PDM_bi.',[],1);

for ii =1:length(PDM_bi_reshape)
    if(PDM_bi_reshape(ii) == 0)
        PDM_bi_reshape(ii) = -1;
    end
end

% F_PDMraw = abs(fft(PDM_bi_reshape));
% figure()
% plot(F_PDMraw)

fs = 8000;
t = ( 0:length(PCM)-1 )/8000;
figure(111)
plot(t,PCM/400)

% F_PCM = abs(fft(PCM/400));
% figure()
% plot(F_PCM)

% hap = dsp.AudioPlayer('SampleRate',fs);
% step(hap,PCM/400);


%
%......................... CIC filter decimation
% R = 64;          %коэф децим
% D = 1;          %задержка
% N = 5;          %порядок фильтра
% 
% cicdec = dsp.CICDecimator(R,D,N);  
% % freqz(cicdec);
% out_CIC = cicdec(PDM_bi_reshape);


% F_out_CIC = abs(fft(out_CIC));
% figure(666)
% plot(F_out_CIC)
%......................... ИЛИ LP filter
out_FIR = decimate(PDM_bi_reshape,64,256,'fir');
%........................HP filter

HP = fir1(64,300/8000,'high');
outCIC_HP = filter(HP,1,out_FIR);
figure(222)
plot(outCIC_HP)

% F_outCIC_HP = abs(fft(outCIC_HP));
% figure(777)
% plot(F_outCIC_HP)

hap = dsp.AudioPlayer('SampleRate',fs);
%%
step(hap,PCM/400);
%%
% step(hap,outCIC_HP/1e7);
step(hap,outCIC_HP*1e2);
%......................... FIR decimation
% 
% outCIC_FIR1 = decimate(out_CIC,2,30,'fir');
% outCIC_FIR2 = decimate(outCIC_FIR1,2,30,'fir');
% outCIC_FIR3 = decimate(outCIC_FIR2,2,30,'fir');
% 
% F_outCIC_FIR2 = abs(fft(outCIC_FIR3));
% % figure()
% % plot(F_outCIC_FIR2)
% %......................... Low pass filter
% LPez = fir1(10,0.2/48,'low');
% 
% 
% figure();
% freqz(LPez,1,512)
% outCIC_FIR3__LP = filter(LPez,1,outCIC_FIR3);
% 
% 
% F_outCIC_FIR3__LP = abs(fft(outCIC_FIR3__LP));
% F_outCIC_FIR3__LP(1) = 0;
% figure()
% plot(F_outCIC_FIR3__LP)
% %
% LP = fir1(16,1/72,'low');
% figure();
% freqz(LP,1,512)
% HP = fir1(255,200/8000,'high');
% 
% 
% out_LP1 = 64*filter(LP,1,PDM_bi_reshape);
% 
% F_out_LP1 = abs(fft(out_LP1));
% F_in = abs(fft(PDM_bi_reshape));
% 
% figure(11);
% plot(F_in)
% 
% figure(12);
% plot(F_out_LP1)
% 
% 
% 
% % out_LP2 = 4*filter(LP,1,out_LP1(1:4:end));
% % out_LP3 = 4*filter(LP,1,out_LP2(1:4:end));
% 
% out_LP3_HP1 = filter(HP,1,out_LP1(1:64:end));
% 
% 
% my_out_filted = out_LP3_HP1;
% % figure();
% % plot(my_out_filted)
% 
% F_PDM_my_filted = abs(fft(my_out_filted));
% figure(666)
% plot(F_PDM_my_filted)
% 
% 
% 
% 
% 
% % simout = conv(simout,LP);
% % simout(end-(length(LP)-1)/2:end) = [];
% % simout(1:(length(LP)-1)/2) = [];
% 
% % simout = conv(simout,LP);
% % simout(end-(length(LP)-1)/2:end) = [];
% % simout(1:(length(LP)-1)/2) = [];
% % 
% % simout = conv(simout,LP);
% % simout(end-(length(LP)-1)/2:end) = [];
% % simout(1:(length(LP)-1)/2) = [];
% 
% 
% 
