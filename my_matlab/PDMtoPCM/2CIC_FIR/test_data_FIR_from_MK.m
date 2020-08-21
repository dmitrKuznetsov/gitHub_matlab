clc; clear; close all;
load('PCM_raw.mat','PCM')
load('PDM_raw.mat','PDM')
pcm_err_ms = 8*50; % первые 50 мс отрезаем 
pdm_err_ms = 64*50;
PCM(1:pcm_err_ms) = [];
PDM(1:pdm_err_ms) = [];

PDM_bi = fliplr( de2bi(PDM(1:64*20)) );

PDM_reshape = reshape(PDM_bi.',[],1);%(1:512)
y = PDM_reshape;
for ii =1:length(y)
    if(y(ii) == 0)
        y(ii) = -1;
    end
end


out_CIC = CIC_nmy_N4(y,32,32);  
out_CIC = floor(out_CIC/(16));


Str = 10;
fileID = fopen('PDM_dec32data.txt','w');
for ii = 1:length(out_CIC)

    fprintf(fileID,'%d, ',out_CIC(ii));
    if (mod(ii,Str)==0)
        fprintf(fileID,'\n');
    end
end

% load('out_CIC_MK.txt')
% figure(1)
% hold on 
% plot(out_CIC)
% plot(out_CIC_MK)
%%
load('ImResp_CompFIR.mat')
out_CIC_FIR = filter(bhi,1,out_CIC);
out_CIC_FIR = out_CIC_FIR(1:2:end);
load('out_CIC_MK.txt')

figure(1)
hold on 
plot(out_CIC_MK/max(abs(out_CIC_MK)))
plot(out_CIC_FIR/max(abs(out_CIC_FIR)))

% figure(45)
% plot(out_CIC_FIR/max(abs(out_CIC_FIR))-out_CIC_MK.'/max(abs(out_CIC_MK)))













