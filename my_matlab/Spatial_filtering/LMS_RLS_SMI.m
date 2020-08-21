clc;clear;close all;
c = 343;
Nel = 2;

f_required = 2000;
lambda = c/f_required;
d = lambda/2;

fs = 8e3;
t_duration = 2;
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);
 
%.................................Сигнал......................................
dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
    'SamplesPerFrame',NSampPerFrame);
sig = step(dftFileReader);

FFsign = fft(sig);
FFsign(1) = 0;
sig = ifft(FFsign);

sigPwr = mean((sig).^2);

%...............  ................Помеха................................
NoisdftFileReader = dsp.AudioFileReader('nois.wav',...%FemaleSpeech-16-8-mono-3secs.wav   nois.wav   my_nois  'Laughter-16-8-mono-4secs.wav'
    'SamplesPerFrame',NSampPerFrame);
noise = 5*step(NoisdftFileReader);

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

F_sign(end,:) = real(F_sign(end,:));
F_nois(end,:) = real(F_nois(end,:));

for i = 1:Nel
    F_sign(pos_len+1:2*pos_len-2,i) = conj(F_sign(pos_len-1:-1:2,i));    
    F_nois(pos_len+1:2*pos_len-2,i) = conj(F_nois(pos_len-1:-1:2,i));
end

t_sign_in = ifft(F_sign,[]);
t_nois_in = ifft(F_nois,[]);
s_in = t_sign_in + t_nois_in;
%%
len_posled = 50:5:300;
% len_posled = length(s_in(:,1));
% LMS 
L = 32;
K = Nel;
mu = 0.05;

for lp = 1: length(len_posled)
    lp

    %LMS
    [y_LMS,W_LMS] = my_Frost(s_in(1:len_posled(lp),:), L, K, mu);
    %RLS 
    [y_LC_RLS, W_RLS] = func_LC_RLS_broadband(s_in(1:len_posled(lp),:),K,L,1);
    
    %SMI
    W_SMI = SMI_func(s_in(1:len_posled(lp),1),s_in(1:len_posled(lp),2),L);
    
    %Сигнал
    [sign_LMS] = apply_coef(t_sign_in, W_LMS, L, K);
    [sign_RLS] = apply_coef(t_sign_in, W_RLS, L, K);
    [sign_SMI] = apply_coef(t_sign_in, W_SMI, L, K);
    %Помеха
    [nois_LMS] = apply_coef(t_nois_in, W_LMS, L, K);
    [nois_RLS] = apply_coef(t_nois_in, W_RLS, L, K);
    [nois_SMI] = apply_coef(t_nois_in, W_SMI, L, K);
    
    SNR_LMS_db(lp) = pow2db(mean(sign_LMS.^2)/mean(nois_LMS.^2));
    SNR_RLS_db(lp) = pow2db(mean(sign_RLS.^2)/mean(nois_RLS.^2));
    SNR_SMI_db(lp) = pow2db(mean(sign_SMI.^2)/mean(nois_SMI.^2));
end
%%
figure(23241)
% title('2 микро, реалная акустическая помеха')
title('2 микро, речевая помеха')
hold on
plot(len_posled, SNR_LMS_db-SNR_true_db, 'b')
plot(len_posled, SNR_RLS_db-SNR_true_db, 'r')
plot(len_posled, SNR_SMI_db-SNR_true_db, 'k')
legend('LMS','RLS','SMI')
xlabel('длина обучающей выборки, сэмпл')
ylabel('выигрыш в ОСШ, дБ')





