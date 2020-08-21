clc; clear; close all;
load('PCM_raw.mat','PCM')
load('PDM_raw.mat','PDM')
pcm_err_ms = 8*50; % первые 50 мс отрезаем 
pdm_err_ms = 64*50;
PCM(1:pcm_err_ms) = [];
PDM(1:pdm_err_ms) = [];

PDM_bi = fliplr( de2bi(PDM(1:64*20)) );

write_buff = reshape(PDM_bi.',[],1);

Str = 10;
fileID = fopen('PDMdata.txt','w');
for ii = 1:length(write_buff)

    fprintf(fileID,'0x%x, ',write_buff(ii));
    if (mod(ii,Str)==0)
        fprintf(fileID,'\n');
    end
end

%%
fclose(fileID);

% iir filter
[b,a] = butter(2,0.5); 

b = round(32767*b)/32767;
a = round(32767*a)/32767;


biquad = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',[b a],'ScaleValues',1);
% freqz(biquad)
out_biquad = biquad(write_buff);
% plotSpectr(out_biquad,9999,8000*64)
%% Фильтрация 
load('OutData.txt')
% out = Hd(write_buff);
% plotSpectr(out,11,8000*64)


% figure(1)
% hold on
% plot(out)
% plot(VarName1)

plotSpectr(write_buff,33,8000*64)
hold on
plotSpectr(OutData,33,8000*64)%out_biquad   OutData

f = (0:length(sig)-1)/length(sig)*fs/2;
Spec = abs(fft(sig));
figure(numFig)
plot(f,20*log10(Spec/max(Spec)) )




% LP = fir1((256-1),300/8000,'low');
% LP = round(LP*1e5);
% outCIC_HP = filter(LP,1,write_buff);
% figure(222)
% plot(outCIC_HP)
% 
% %%
% fileID = fopen('LPdata.txt','w');
% for ii = 1:length(LP)
% 
%     fprintf(fileID,'%d, ',LP(ii));
%     if (mod(ii,Str)==0)
%         fprintf(fileID,'\n');
%     end
% end
% fclose(fileID);