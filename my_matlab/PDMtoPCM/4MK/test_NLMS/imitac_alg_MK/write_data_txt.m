% clc; clear; close all;
write_buff = mic2(1:8000);
Str = 10;
fileID = fopen('voic2_f32.txt','w');
for ii = 1:length(write_buff)

    fprintf(fileID,'%E, ',write_buff(ii));
    if (mod(ii,Str)==0)
        fprintf(fileID,'\n');
    end
end
fclose(fileID);