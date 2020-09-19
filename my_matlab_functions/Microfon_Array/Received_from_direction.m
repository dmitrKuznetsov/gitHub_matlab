function t_sign_in = Received_from_direction(signal, Nel, teta, d, fs)
if(~iscolumn(signal))
    signal = signal.';
    if(~iscolumn(signal))
        warning('В функцию подана матрица');
    end
end
c = 343;
    
f = (0:length(signal)-1)/length(signal) * fs;

f_pos = f(2:length(f)/2+1);                                       
num_bin = length(f_pos);                                                 
freq_per_bin = max(f_pos)/num_bin;                                
sample_per_bin = length(f_pos)/num_bin;                           

pos_len = length(f_pos)+1;
FFsign = fft(signal);
F_sign = zeros(pos_len,Nel);                                       

for i_bin = 1:num_bin
    f_center = (2*i_bin-1)*freq_per_bin/2;                                  % (f_min+f_max)/2                       
    index = (i_bin-1)*sample_per_bin+2  :  i_bin*sample_per_bin+1;
    for i = 1:Nel
        F_sign(index, i) = FFsign(index)*exp( -1i*2*pi*f_center/c*d*(i-1)*sind(teta) );
    end
end
F_sign(end) = real(F_sign(end));

for i = 1:Nel
    F_sign(pos_len+1:2*pos_len-2,i) = conj(F_sign(pos_len-1:-1:2,i));    
end
t_sign_in = ifft(F_sign);

if(~isreal(t_sign_in))
    warning('Выход комплексный');
end

end

