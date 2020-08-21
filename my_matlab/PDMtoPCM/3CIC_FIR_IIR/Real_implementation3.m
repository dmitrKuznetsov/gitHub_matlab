% PCM = VarName1;
% PDM = VarName2;
% save('PCM_raw_2.mat','PCM')
% save('PDM_raw_2.mat','PDM')
%%
clc; clear; close all;
load('PCM_raw_2.mat','PCM')
load('PDM_raw_2.mat','PDM')
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

R = 32;          %коэф децим
D = 2;          %задержка
N = 4;          %порядок фильтра
% .........................nmy CIC filter decimation......................

out_CIC = floor(CIC_nmy_N4_v2(y,R,32)/8);  
% plotSpectr(out_CIC_test,666,8000*64/R)

%% .........................correctir filtr......................
in = [0 1 ones(1,16*1200-2)];                   %единичное возмущение
Imp_resp = CIC_nmy_N4_v2(in,R,32);            

Freq_resp = abs(fft(Imp_resp));

Freq_resp(1) = 0;
Freq_resp = Freq_resp/max(Freq_resp);
Freq_resp(1) = 1;

Fresp_ideal = [ones(1,length(Freq_resp)/4 ) zeros(1,length(Freq_resp)/2 ) ones(1, length(Freq_resp)/4)];


Fresp_correct = Fresp_ideal./Freq_resp;
Fresp_correct(length(Fresp_correct)/2+1:end) = [];
% figure(777)
% plot((Fresp_correct))

f = (0:length(Fresp_correct))/length(Fresp_correct);

bhi = fir2(32-1,f,[Fresp_correct 0] );                %!!!!!
bhi =  round(bhi*10000);
% save('ImResp_CompFIR.mat', 'bhi')
out_CIC_corFIR = filter(bhi,1,out_CIC);


% LP = fir1(32-1,1/2, 'low');
% out_CIC_corFIR = filter(LP,1,out_CIC);

% plotSpectr(out_correct_FIR,789,8000*2)

out_CIC_corFIR_dec = out_CIC_corFIR(1:2:end);
% plotSpectr(out_correct_FIR_decim,987,8000)
%% HP filter 
% HP = fir1(256,150/8000,'high');
% freqz(HP)
% out_CIC_corFIR_HP = filter(HP,1,out_CIC_corFIR_dec );

[b,a] = butter(2,150/8000,'high'); 
b = round(32767*b)/32767;
a = round(32767*a)/32767;

biquad = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',[b a],'ScaleValues',1);
out_CIC_corFIR_HP = biquad(out_CIC_corFIR_dec.');



plotSpectr(out_CIC_corFIR_HP,1000,8000)
hold on
plotSpectr(PCM,1000,8000)
%%
figure(167)
hold on 
plot(PCM/max(abs(PCM)))
plot(out_CIC_corFIR_HP/max(abs(out_CIC_corFIR_HP)))


hap = dsp.AudioPlayer('SampleRate',fs);
% %% Слушаем обычный
% load('PCM_normal.mat','PCM_normal')
% step(hap,PCM_normal/max(PCM_normal));
% %% Слушаем скомпенсированный
% step(hap,out_CIC_corFIR_HP/max(out_CIC_corFIR_HP));
% %%
% 
% step(hap,PCM/max(PCM));






