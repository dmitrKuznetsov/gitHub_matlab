clc;clear;close all;
load('IIR_4order_2secBiquad_120Fc.mat')
% load('IIR_6order_3secBiquad_120Fc.mat')
biquad1 = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',SOS,'ScaleValues',1);

freqz(biquad1)

SOS2 = round(SOS/4*32767);
biquad2 = dsp.BiquadFilter('Structure','Direct form I',...
    'SOSMatrix',SOS2,'ScaleValues',1);

freqz(biquad2)







