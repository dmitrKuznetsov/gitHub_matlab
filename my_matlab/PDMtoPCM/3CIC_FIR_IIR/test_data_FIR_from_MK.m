clc; clear; close all;
load('PCM_raw.mat','PCM')
load('PDM_raw.mat','PDM')
pcm_err_ms = 8*50;                  % первые 50 мс отсекаем 
pdm_err_ms = 64*50;
PCM(1:pcm_err_ms) = [];
PDM(1:pdm_err_ms) = [];
Write_data(PDM(1:64*20),'PDMraw1.txt');% 20мс записываем
Write_data(PDM(64*20+1:64*40),'PDMraw1.txt');% 20мс записываем
PDM_bi = fliplr( de2bi(PDM(1:64*20)) );

PDM_reshape = reshape(PDM_bi.',[],1);%(1:512)
y = PDM_reshape;
for ii =1:length(y)
    if(y(ii) == 0)
        y(ii) = -1;
    end
end


out_CIC = CIC_nmy_N4_v2(y,32,32);  
out_CIC = floor(out_CIC/(16));




% load('out_CIC_MK.txt')
% figure(1)
% hold on 
% plot(out_CIC)
% plot(out_CIC_MK)
%% FIR 
load('ImResp_CompFIR.mat')
out_CIC_FIR = filter(bhi,1,out_CIC);
out_CIC_FIR = out_CIC_FIR(1:2:end);
out_CIC_FIR = floor(out_CIC_FIR/(2^15));
%% IIR
Str = 10;
fileID = fopen('out_CIC_FIR.txt','w');
for ii = 1:length(out_CIC_FIR)

    fprintf(fileID,'%d, ',out_CIC_FIR(ii));
    if (mod(ii,Str)==0)
        fprintf(fileID,'\n');
    end
end


%% IIR biquad
% HP = fir1(256-1,150/8000,'high');
% out_CIC_FIR_IIR = filter(HP,1,out_CIC_FIR.');

[b,a] = butter(2,150/8000,'high'); 
SOS = [b a];
% SOS = round(32767*SOS)/32767;
% SOS = [1029936539	-2059873078 	1029936539   (2^31-1)/2 	-2058085216		987919116]*2/(2^31-1);
biquad = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',SOS,'ScaleValues',1);
out_CIC_FIR_IIR = biquad(out_CIC_FIR.');


% load('IIR_4order_2secBiquad_120Fc.mat')
% load('IIR_6order_3secBiquad_120Fc.mat')
% biquad1 = dsp.BiquadFilter('Structure','Direct form I',...
%     'SOSMatrix',SOS,'ScaleValues',1);
% out_CIC_FIR_IIR = biquad1(out_CIC_FIR.');



load('out_CIC_FIR_IIR_MK.txt')

figure(1)
hold on 
plot(out_CIC_FIR_IIR_MK)
plot(out_CIC_FIR_IIR)

%%
plotSpectr(out_CIC_FIR_IIR,1000,8000)
hold on
plotSpectr(out_CIC_FIR_IIR_MK,1000,8000)


%%










