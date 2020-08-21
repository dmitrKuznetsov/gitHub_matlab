clc;clear;close all;
c = 343;
Nel = 8;

f_required = 2000;
lambda = c/f_required;
d = lambda/2;

fs = 8e3;
t_duration = 2;
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);
 
% hap = dsp.AudioPlayer('SampleRate',fs);
%.................................Сигнал......................................
dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
    'SamplesPerFrame',NSampPerFrame);
sig = 0.05*step(dftFileReader);

FFsign = fft(sig);
FFsign(1) = 0;
sig = ifft(FFsign);

sigPwr = mean((sig).^2);

%...............  ................Помеха................................

fsin = 1000;
noise = sin(2*pi*fsin*t);

FFnois = fft(noise);
FFnois(1) = 0;
noise = ifft(FFnois);

noisPwr = mean((noise).^2);

SNR_true = sigPwr/noisPwr;
SNR_true_db = pow2db(SNR_true);

tetaSign = 0;                    %[phi;teta]
tetaNois = 55;
% ................................Прием................................
f = (0:length(sig(:,1))-1)/t_duration;
% Обрабатываю положительную часть спектра
f_pos = f(2:length(f)/2+1);                                       

num_bin_in = length(f_pos);                                                   
freq_per_bin = max(f_pos)/num_bin_in;                                
sample_per_bin = length(f_pos)/num_bin_in;                           

pos_len = length(f_pos)+1;
F_sign = zeros(pos_len,Nel);                                       
F_nois = zeros(pos_len,Nel);

for i_bin = 1:num_bin_in
    f_00 =(2*i_bin-1)*freq_per_bin/2;                             
    index = (i_bin-1)*sample_per_bin+2:i_bin*sample_per_bin+1;
    for i=1:Nel
        F_sign(index,i) = FFsign(index)*exp( 1i*2*pi*f_00/c*d*(i-1)*sind(tetaSign) );
        F_nois(index,i) = FFnois(index)*exp( 1i*2*pi*f_00/c*d*(i-1)*sind(tetaNois) );
    end
end

for i = 1:Nel
F_sign(pos_len+1:2*pos_len-2,i) = conj(F_sign(pos_len-1:-1:2,i));    
F_nois(pos_len+1:2*pos_len-2,i) = conj(F_nois(pos_len-1:-1:2,i));
end

t_sign_in = ifft(F_sign,[],1,'symmetric');
t_nois_in = ifft(F_nois,[],1,'symmetric');
s_in = t_sign_in + t_nois_in;
% ................................Обработка.............................
J = 25;
K = Nel;
mu = 0.005;
[y1,W1] = my_Frost(s_in, J, K, mu);

%%
W = zeros(K,J);
X = zeros(K,J);
y2 = zeros(length(s_in),1);
if(J == 1)
    C = 1;
else
    C = fir1(J-1,0.99,'low',chebwin(J,30)).';
end

for ii = 1:length(W(:,1))
    W(ii,:) = 1/K*C.';
end

for iter = 1:length(s_in)
    X(:,J) = [];
    X = [s_in(iter, :).' X];
    y2(iter) = sum(sum(W.*X));
    W = W - mu*y2(iter)*X;
    W = W + 1/K*(C.'-ones(1,K)*W);

end

figure(3)
hold on
plot(t,y1-y2);


