% PCM = VarName1;
% PDM = VarName2;
% save('PCM_raw_2.mat','PCM')
% save('PDM_raw_2.mat','PDM')
%%
clc;clear;close all;
% load('PDM_MK.txt')
% Write_data(PDM_MK,'PDMraw123.txt');
% load('DATA\PCM_raw.mat','PCM')
% load('DATA\PDM_raw.mat','PDM')

load('DATA\3PCM.log')
load('DATA\3PDM.log')
PDM = X3PDM;
PCM = X3PCM;

PDM_bi = fliplr( de2bi(PDM) );

PDM_reshape = reshape(PDM_bi.',[],1);
for ii =1:length(PDM_reshape)
    if(PDM_reshape(ii) == 0)
        PDM_reshape(ii) = -1;
    end
end

PCM_matlab = PDM2PCM(PDM_reshape);
% load('PCM_MK.txt')


figure(1)
hold on

plot(PCM_matlab/max( abs(PCM_matlab(300:end)) ))    %   /max( abs(PCM_matlab(300:end)) )
plot(PCM/max(abs(PCM(300:end))))
legend('matlab','MK')
ylim([-2 2])





plotSpectr(PCM_matlab,1000,8000)
hold on
plotSpectr(PCM,1000,8000)
%%
% sound(PCM_matlab/max(abs(PCM_matlab)),8000) 




