clc;clear;close all;
load('noise_voice34mic_4.mat','VarName1');


X = VarName1;
voic1 = X(1:length(X)/2);
voic2 = X(length(X)/2+1:end);
voic1(1:400) = [];
voic2(1:400) = [];

voic1(16000+1:end) = [];
voic2(16000+1:end) = [];
%%
voic1 = voic1+10000;
voic2 = voic2+10000;


intelhex_gen_my_v1(voic1);
%%
intelhex_gen_my_v2(voic2);