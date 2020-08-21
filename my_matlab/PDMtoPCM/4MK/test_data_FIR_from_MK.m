clc; clear; close all;
load('PDM_MK.txt')
load('PCM_MK.txt')
% PDM_MK = (2^8-1)*ones(640,1);

% Write_data(PDM_MK,'PDMraw1.txt');

PDM_bi = fliplr( de2bi(PDM_MK) );

PDM_reshape = reshape(PDM_bi.',[],1);%(1:512)
y = PDM_reshape;
for ii =1:length(y)
    if(y(ii) == 0)
        y(ii) = -1;
    end
end


out_CIC = CIC_nmy_N4_v2(y,32,32);  
% load('outCIC_MK.txt')
% 
% plot(666)
% hold on
% plot(out_CIC)
% plot(outCIC_MK)

%% FIR 
load('ImResp_CompFIR.mat')
out_CIC_FIR = filter(bhi,1,out_CIC);
out_CIC_FIR = out_CIC_FIR(1:2:end);
out_CIC_FIR = floor(out_CIC_FIR/(2^31));

load('out_CIC_FIR_MK.txt')

% figure(777)
% hold on
% plot(out_CIC_FIR)
% plot(out_CIC_FIR_MK)



%% IIR biquad


[b,a] = butter(2,150/8000,'high'); 
SOS = [b a];
SOS = round((2^31-1)*SOS)/(2^31-1);

biquad = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',SOS,'ScaleValues',1);
PCM = biquad(out_CIC_FIR.');




figure(1)
hold on 
plot(PCM_MK)
plot(PCM)

% %%
plotSpectr(PCM_MK,1000,8000)
hold on
plotSpectr(PCM,1000,8000)











