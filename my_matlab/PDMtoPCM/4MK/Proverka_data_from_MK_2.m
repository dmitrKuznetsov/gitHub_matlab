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

% load('DATA\3PCM.log')
% load('DATA\3PDM.log')
% PDM = X3PDM;
% PCM = X3PCM;


load('PDM.txt');                        %DATA\
PDM(1:mod(length(PDM),64)) = [];
% PDM = X1PDM;
PDM_bi = fliplr( de2bi(PDM) );

PDM_reshape = reshape(PDM_bi.',[],1);
for ii =1:length(PDM_reshape)
    if(PDM_reshape(ii) == 0)
        PDM_reshape(ii) = -1;
    end
end

PCM_my = PDM2PCM(PDM_reshape);
PCM_matlab = PDM2PCM_matlab(PDM_reshape);
% load('PCM_MK.txt')


figure(1)
hold on
plot(PCM_my/max(PCM_my))
plot(PCM_matlab/max(PCM_matlab))                                                                                                                                                                                                                                                     
legend('MK','matlab')

plotSpectr(PCM_my,1000,8000)
hold on
plotSpectr(PCM_matlab,1000,8000)

%%
sound(PCM_matlab/max(PCM_matlab),8000) 




