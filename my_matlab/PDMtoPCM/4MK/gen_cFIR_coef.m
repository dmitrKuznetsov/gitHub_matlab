clc;clear;close all;


in = [0 1 ones(1,16*1200-2)];                   %единичное возмущение
Imp_resp = CIC_nmy_N4_v2(in,32,32);            

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
bhi =  round(bhi*(2^31-1));
figure(888)
freqz(bhi)

Write_data(bhi,'coef_cFIR.txt')


% save('ImResp_CompFIR.mat', 'bhi')












