function [PCM] = PDM_uint8_2PCM_my(PDM, str, numMic)


if(str == 'lr')
    PDM_bi = fliplr( de2bi(PDM) );
else
    PDM_bi = ( de2bi(PDM) );
end
PDM_reshape = reshape(PDM_bi.',[],1);

if(numMic == 1)
    PDM_reshape_1mic = PDM_reshape;
else
    PDM_reshape_1mic = PDM_reshape(1:2:end);
end

for ii = 1:length(PDM_reshape_1mic)
    if(PDM_reshape_1mic(ii) == 0)
        PDM_reshape_1mic(ii) = -1;
    end
end
% for ii = 1:length(PDM_reshape_1mic)
%     if(PDM_reshape_1mic(ii) == 0)
%         PDM_reshape_1mic(ii) = -128;
%     else
%         PDM_reshape_1mic(ii) = 127;
%     end
% end

PDM_reshape_1mic(1:mod(length(PDM_reshape_1mic),32)) = [];

R = 32;         %коэф децим
out_CIC = CIC_nmy_N4_v2(PDM_reshape_1mic,R,32);  
out_CIC = floor(out_CIC/(16));

%% .........................correctir filtr......................
load('ImResp_CompFIR.mat');
out_CIC_FIR = filter(bhi,1,out_CIC);
out_CIC_FIR = out_CIC_FIR(1:2:end);
out_CIC_FIR = floor(out_CIC_FIR/(2^15));
%% HP filter 

[b,a] = butter(2,150/8000,'high'); 
SOS = [b a];
SOS = round(32767*SOS)/32767;
biquad = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',SOS,'ScaleValues',1);
PCM = biquad(out_CIC_FIR.');
% HP = fir1(256,300/8000,'high');
% PCM = filter(HP,1,out_CIC_FIR);

end

