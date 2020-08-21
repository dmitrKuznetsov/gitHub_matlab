function [PCM] = PDM2PCM_matlab(PDM)

R = 32;          %коэф децим
D = 2;          %задержка
N = 4;          %порядок фильтра

cicdec1 = dsp.CICDecimator(R,D,N);

out_CIC = cicdec1(PDM);  
out_CIC = floor(out_CIC/(16));

%% .........................correctir filtr......................
load('ImResp_CompFIR.mat');
out_CIC_FIR = filter(bhi,1,out_CIC);
out_CIC_FIR = out_CIC_FIR(1:2:end);
out_CIC_FIR = floor(out_CIC_FIR/(2^15));
%% HP filter 

% [b,a] = butter(2,150/8000,'high'); 
% SOS = [b a];
% SOS = round(32767*SOS)/32767;
% biquad = dsp.BiquadFilter('Structure','Direct form I',...
%     'SOSMatrix',SOS,'ScaleValues',1);
% PCM = biquad(out_CIC_FIR.');
HP = fir1(256,300/8000,'high');
PCM = filter(HP,1,out_CIC_FIR);

end

