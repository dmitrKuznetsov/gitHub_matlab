clc;clear;close all;
load('PDM_MK.txt')
% Write_data(PDM_MK,'PDMraw123.txt');

PDM_bi = fliplr( de2bi(PDM_MK) );

PDM_reshape = reshape(PDM_bi.',[],1);
for ii =1:length(PDM_reshape)
    if(PDM_reshape(ii) == 0)
        PDM_reshape(ii) = -1;
    end
end

PCM_matlab = PDM2PCM(PDM_reshape);
%load('PCM_MK.txt')


figure(1)
hold on
%%plot(PCM_MK/max(PCM_MK))
plot(PCM_matlab/max(PCM_matlab))


%plotSpectr(PCM_MK,1000,8000)
hold on
plotSpectr(PCM_matlab,1000,8000)






