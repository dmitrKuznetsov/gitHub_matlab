clc; clear; close all;
load('PCM_raw.mat','PCM')
load('PDM_raw.mat','PDM')
pcm_err_ms = 8*50; % первые 50 мс отрезаем 
pdm_err_ms = 64*50;
PCM(1:pcm_err_ms) = [];
PDM(1:pdm_err_ms) = [];

PDM_bi = fliplr( de2bi(PDM) );

PDM_reshape = reshape(PDM_bi.',[],1);%(1:512)
y = PDM_reshape;
for ii =1:length(y)
    if(y(ii) == 0)
        y(ii) = -1;
    end
end



fs = 8000;
PCM = PCM/max(PCM);
t = ( 0:length(PCM)-1 )/8000;

R = 32;         %коэф децим
D = 2;          %задержка
N = 4;          %порядок фильтра

out_CIC = CIC_nmy_N4(y,R,38);

%% LP
dec = 2;
LP = fir1(32-1,1/dec, 'low');
% freqz(LP,1,512)

out_CIC_LP = filter(LP,1,out_CIC);
out_CIC_LP = out_CIC_LP(1:dec:end);
%% HP
HP = fir1(256-1,150/8000,'high');
% freqz(HP)
out_CIC_LP_HP = filter(HP,1,out_CIC_LP );


figure(1)
plot(out_CIC_LP_HP)



plotSpectr(PCM,1000,8000)
hold on
plotSpectr(out_CIC_LP_HP,1000,8000)
%% Слушаем ....
hap = dsp.AudioPlayer('SampleRate',fs);
step(hap,out_CIC_LP_HP.'/max(out_CIC_LP_HP));
%% Корректированный
load('voic_corFIR.mat','voic_corFIR')
step(hap,voic_corFIR);
%%
step(hap,PCM/max(PCM));




% plotSpectr(out_CIC_LP_HP,1001,8000)
% hold on
% plotSpectr(voic_corFIR,1001,8000)