clc; clear; close all;
load('PCM_raw.mat','PCM')
load('PDM_raw.mat','PDM')
pcm_err_ms = 8*50; % первые 50 мс отрезаем 
pdm_err_ms = 64*50;
PCM(1:pcm_err_ms) = [];
PDM(1:pdm_err_ms) = [];

PDM_bi = fliplr( de2bi(PDM) );

y = reshape(PDM_bi.',[],1);%(1:512)

for ii =1:length(y)
    if(y(ii) == 0)
        y(ii) = -1;
    end
end



fs = 8000;
PCM = PCM/max(PCM);
t = ( 0:length(PCM)-1 )/8000;

R = 32;          %коэф децим
D = 2;          %задержка
N = 4;          %порядок фильтра
cicdec1 = dsp.CICDecimator(R,D,N);  


% .........................nmy CIC filter decimation......................
out_CIC_test = CIC_nmy_N4(y,R,32);   %28 34 37

out_CIC_mat = cicdec1(y);                                   


plotSpectr((out_CIC_mat),666,8000*64/R)
hold on
plotSpectr(out_CIC_test,666,8000*64/R)

%%
figure(777);
hold on
plot(out_CIC_mat)
plot(out_CIC_test)





