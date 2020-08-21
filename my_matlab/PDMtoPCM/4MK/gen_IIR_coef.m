clc;clear;close all;
[b,a] = butter(2,150/8000,'high'); 
SOS = [b a];
SOS = SOS/2;
SOS = round((2^31-1)*SOS);

load('ImResp_CompFIR.mat')











