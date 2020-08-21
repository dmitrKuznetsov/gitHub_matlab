clc;clear; 
fs = 8e3;
t_duration = 2;
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);
 

dftFileReader = dsp.AudioFileReader('nois.wav',...
    'SamplesPerFrame',NSampPerFrame);
nois = step(dftFileReader);
nois_mass = []
for i = 1: 60
    nois_mass = [nois_mass; nois];
    
end
audiowrite('nois_mass.wav',nois_mass,8e3)



