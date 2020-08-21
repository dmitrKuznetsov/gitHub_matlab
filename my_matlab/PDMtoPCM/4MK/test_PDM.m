% PCM = VarName1;
% PDM = VarName2;
% save('PCM_raw_2.mat','PCM')
% save('PDM_raw_2.mat','PDM')
%%
% clc;clear;
% 
% load('DATA\3PDM.log')
% PDM = X3PDM;

% load('DATA\3PCM.log')
% PCM = X3PCM;





% PDM = Value;
PDM = VarName1(1:length(VarName1)/2);

PCM_my = PDM_uint8_2PCM_my(PDM, 'lr', 2);
PCM_matlab = PDM_uint8_2PCM_matlab(PDM, 'lr', 2);


%%
close all;
figure(1)
hold on
plot(PCM_my/max(PCM_my)/2)
plot(PCM_matlab/max(PCM_matlab))                                                                                                                                                                                                                                                     
legend('MK','matlab')
ylim([-1 1])
% plotSpectr(PCM_my,1000,8000)
% hold on
% plotSpectr(PCM_matlab,1000,8000)

%%
%  sound(PCM_matlab/max(PCM_matlab),8000) 
 sound(PCM_my/max(PCM_my),8000) 



