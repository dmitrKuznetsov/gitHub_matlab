clc;clear;close all;
load('ImResp_CompFIR.mat', 'bhi')
fileID = fopen('ImResp_CompFIR.txt','w');
Str = 10;
for ii = 1:length(bhi)
    
    fprintf(fileID,'%d, ',bhi(ii));
    
    
    
    if (mod(ii,Str)==0)
        fprintf(fileID,'\n');
    end
end