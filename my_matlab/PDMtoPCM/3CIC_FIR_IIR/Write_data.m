function Write_data(write_buff,name)

str_num_el = 10;
fileID = fopen(name,'w');
for ii = 1:length(write_buff)
    fprintf(fileID,'0x%x, ',write_buff(ii));
%     fprintf(fileID,'%d, ',write_buff(ii));
    if (mod(ii,str_num_el)==0)
        fprintf(fileID,'\n');
    end
end
end

