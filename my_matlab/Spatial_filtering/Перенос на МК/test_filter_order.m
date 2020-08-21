clc;clear;close all;
c = 343;
Nel = 2;

f_required = 2000;
d = c/f_required/2;

fs = 8e3;
t_duration = 2;
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);
 
% hap = dsp.AudioPlayer('SampleRate',fs);
%.................................Сигнал......................................
dftFileReader = dsp.AudioFileReader('dft_voice_8kHz.wav',...
    'SamplesPerFrame',NSampPerFrame);
sig = 400*step(dftFileReader);

FFsign = fft(sig);
FFsign(1) = 0;
sig = ifft(FFsign);

sigPwr = mean((sig).^2);

%...............  ................Помеха................................
% SNR_true=0.01;
% SNR_true_db=pow2db(SNR_true);
% noisePwr = sigPwr/SNR_true; 
% noise = sqrt(noisePwr)*randn(length(sig),1); %pinknoise(length(sig), 1)%randn(length(sig),1)
% noise = noise-mean(noise);

NoisdftFileReader = dsp.AudioFileReader('nois.wav',...%FemaleSpeech-16-8-mono-3secs.wav   nois.wav
    'SamplesPerFrame',NSampPerFrame);
noise = 400*step(NoisdftFileReader);

FFnois = fft(noise);
FFnois(1) = 0;
noise = ifft(FFnois);

noisPwr = mean((noise).^2);

SNR_true = sigPwr/noisPwr;
SNR_true_db = pow2db(SNR_true);

tetaSign = 0;                    %[phi;teta]
tetaNois = 40;
% ................................Прием................................
f = (0:length(sig(:,1))-1)/t_duration;
% Обрабатываю положительную часть спектра
f_pos = f(2:length(f)/2+1);                                       

num_bin_in = 1000;                                                   
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

%%
NQWER = 8;
K_ys = 1024;
K = Nel;
% mu = 0.0000005;
% mu = 2000000;
mu = 2^21;

Jmas = 35:35;%50;
for i_J = 1:length(Jmas)
J = Jmas(i_J);
W_mic1_lr = zeros(1,J);
W_mic2_lr = zeros(1,J);


y = zeros(1,length(s_in(:,1)));
if(J == 1)
    C = 1;
else
    C = floor(fir1(J-1,0.99,'low',chebwin(J,30)).'*K_ys).';
end  

C_lr = fliplr(C);
W_mic1_lr = floor(fliplr(1/K*C_lr));
W_mic2_lr = floor(fliplr(1/K*C_lr));
X_mic1 = s_in(:,1).';
X_mic2 = s_in(:,2).';

Niter = fix(length(X_mic1)/NQWER);
data_mic1 = zeros(1, (J-1)+NQWER);
data_mic2 = zeros(1, (J-1)+NQWER);



for i_iter = 0:Niter-1
    index = 1+NQWER*i_iter : NQWER+NQWER*i_iter;
    data_mic1(J:J+NQWER-1) = X_mic1(index);
    data_mic2(J:J+NQWER-1) = X_mic2(index);
    for i_elem = 1:NQWER %!!!!!!!
        y(i_elem+NQWER*i_iter) = 0;
        for i_otvod = 0:J-1
            y(i_elem+NQWER*i_iter) = y(i_elem+NQWER*i_iter) + data_mic1(i_elem+i_otvod)*W_mic1_lr(i_otvod+1)...
                                                            + data_mic2(i_elem+i_otvod)*W_mic2_lr(i_otvod+1);
        end
        for i_otvod = 1:J%!!!
            W_mic1_lr(i_otvod) = W_mic1_lr(i_otvod) - fix(y(i_elem+NQWER*i_iter)*data_mic1(i_elem+(i_otvod-1))/mu);
            W_mic2_lr(i_otvod) = W_mic2_lr(i_otvod) - fix(y(i_elem+NQWER*i_iter)*data_mic2(i_elem+(i_otvod-1))/mu);
            
        end
        for i_otvod = 1:J%!!!
            V_corect(i_otvod) = fix((C_lr(i_otvod) - W_mic1_lr(i_otvod) - W_mic2_lr(i_otvod))/K);
            W_mic1_lr(i_otvod) = W_mic1_lr(i_otvod) + V_corect(i_otvod);
            W_mic2_lr(i_otvod) = W_mic2_lr(i_otvod) + V_corect(i_otvod);
        end
    end
    
    %в преамбулу добавим данные из прошлой итерации
    data_mic1(1:J-1) = data_mic1(end-J+2:end);
    data_mic2(1:J-1) = data_mic2(end-J+2:end);
        
end        
y = fix(y/K_ys);

SNR_outFrost(i_J) = pow2db( mean(t_sign_in(:,1).^2)/(mean((y.').^2) - mean(t_sign_in(:,1).^2) ));

end
%%
% figure(111);
% plot(Jmas,SNR_outFrost)

%%
figure(1);
plot(s_in(:,1))
hold on 
plot(y)
plot(t_sign_in(:,1))
% figure(33);
% plot(t_sign_in(:,1))
% hold on 
% plot(y)
% 
% t_sign_in(1:J,:) = [];
% figure(2);
% plot(y)
% hold on 
% plot(t_sign_in(:,1))
% xlim([8740 8990])




