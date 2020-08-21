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

% for ii =1:length(PDM_bi_reshape)
%     if(PDM_bi_reshape(ii) == 0)
%         PDM_bi_reshape(ii) = -1;
%     end
% end
% plotSpectr(PDM_bi_reshape,000,64*8000)

fs = 8000;
PCM = PCM/max(PCM);
t = ( 0:length(PCM)-1 )/8000;
figure(111)
plot(t,PCM)
% plotSpectr(PCM,222,8000)




HP = fir1(64,300/8000,'high');
%% ......................... CIC filter decimation......................
R = 64;         %коэф децим
D = 1;          %задержка
N = 5;          %порядок фильтра

cicdec = dsp.CICDecimator(R,D,N);  
% freqz(cicdec);
out_CIC = cicdec(PDM_bi_reshape);
% out_CIC(1:300) = [];
%...........................HP filter
outCIC_HP = filter(HP,1,out_CIC);
outCIC_HP(1:300) = [];
outCIC_HP = outCIC_HP/max(outCIC_HP);
figure(333)
plot(outCIC_HP)
% plotSpectr(outCIC_HP,444,8000)

%% ..............................LP filter................................
out_FIR = decimate(PDM_bi_reshape,64,64,'fir');
%...........................HP filter
outFIR_HP = filter(HP,1,out_FIR);
outFIR_HP(1:300) = [];
outFIR_HP = outFIR_HP/max(outFIR_HP);
figure(555)
plot(outFIR_HP)
% plotSpectr(outFIR_HP,666,8000)


% hap = dsp.AudioPlayer('SampleRate',fs);
% %%
% step(hap,PCM);
% %%   CIC
% step(hap,outCIC_HP);
% %%   FIR
% step(hap,outFIR_HP);



