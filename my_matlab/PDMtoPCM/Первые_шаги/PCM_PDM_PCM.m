clc;clear; close all;
fs=8e3;
t_duration=3;
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame=length(t);

dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
    'SamplesPerFrame',NSampPerFrame);
sig = step(dftFileReader);

player = audioplayer(2*sig, 8000);
play(player);
%%
DataSign = [t.', sig];

%%
sim('pcm_pdm_pcm.slx')
%%
k = 16; 
PDM(1:640) = [];

for i = 1:length(PDM)
    if(PDM(i) == -1)
       PDM(i) = 0;
    end
end

PDM_reshape = reshape(PDM(1:floor(length(PDM)/k)*k),[k, floor(length(PDM)/k)]);
PDM_reshape = PDM_reshape.';
PDM_reshape_de = bi2de(PDM_reshape);

data = PDM_reshape_de;
%%













intelhex_gen_my(data(1:32768));

%plot(PDM_reshape_de)