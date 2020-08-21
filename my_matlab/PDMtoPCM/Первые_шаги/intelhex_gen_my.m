%% Intel-hex generation from vector variable
%%%%%%%%%%%%%%%%%%%%%%%
%  by Daniil Smirnov  %
%  updated 2016/08/26 %
%%%%%%%%%%%%%%%%%%%%%%%
% version 1.3
% - added support any rate of data including negative integers with length
% no more than 65535 values
% - warnings and errors added to prevent incorrect hex-file
function [ ] = intelhex_gen_my(data)

if length(data) > 655535
    error('This script does not support vectors with length more than 655535');
end
if rem(data, floor(data)) ~= 0
    error('This script supports only integers.');
end
if isvector(data) == 0
    error('please, use only vector to convert')
end
%prelocation of buffer variables
a = [];
twos = [];
result = [];

fileName = input('Enter valid name for Intel-hex File: ', 's');
if isempty(fileName)
    fileName = 'hex_file';
end
fileName = [fileName '.hex'];
if (exist(fileName, 'file'))
    warning('File with the same name already exist. It will be replaced.')
end

start_string = [':' dec2hex(ceil((length(dec2bin(ceil(max(abs(data))))))/8),2)];
pow = 2^(ceil((length(dec2bin(max(abs(data))))+1)/8)*8);
for ii = 1:length(data)
    
    data_i = dec2hex(uint64(data(ii)),ceil(length(dec2bin(ceil(max(data))))/4));
    if(data(ii)<0)
        data(ii) = data(ii) + pow;
    end
%       +-------------+--------------+---------+----------------------+
%       |    start    |   address    |data type|         data         |
    a = [start_string dec2hex((ii-1)*2+800,4)  '00'  data_i(3:4) data_i(1:2)];% char(10)];
    % writing in text format including CR and LF bytes
    ch = 0;
    % counting checksum:
    for oo = 1:(length(a)-1)/2
        ch = ch + hex2dec(a(2*oo:2*oo+1)); 
    end
    sh = dec2bin(ch,8);
    kk = 1;
    for jj = length(sh) - 7:length(sh)
        if (sh(jj) =='1')
            twos(kk) = '0';
        else
            twos(kk) = '1';
        end
        kk = kk + 1;
    end
    s = dec2hex((bin2dec(char(twos)) + 1),2);
    result = [result a s(end-1:end) char(13) char(10)];
end
% buffer with all result
result = [result ':00000001FF' char(13) char(10)];
% writing all data in file
fileID = fopen(fileName, 'w');
fprintf(fileID, '%c', result);
fclose(fileID);
disp('File generated succesfully')
fprintf('Word size: %d \n', length(dec2bin(ceil(max(abs(data))))));
end