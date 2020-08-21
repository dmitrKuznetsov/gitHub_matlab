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
% ................................Обработка.............................
J = 5;
K = Nel;
mu = 0.005;
[y1,W1] = my_Frost(s_in, J, K, mu);


X = zeros(K*J,1);
y = zeros(K*J,1);
C = zeros(K*J,J);
W = zeros(K*J,1);
I = diag(ones(1,K*J));
P = zeros(K*J,K*J);
F = zeros(K*J,1);                                    
for j = 1:J
    C(:,j) = [zeros(1,(j-1)*K) ones(1,K) zeros(1,J*K-j*K)].';
end
if(J == 1)
    FF = 1;
else
FF = fir1(J-1,0.99,'low',chebwin(J,30)).';
end
% figure(111);
% freqz(FF,1)
F = C/(C.'*C)*FF;                           %!!!C* inv(C.'*C)*FF;
P = I - C/(C.'*C)*C.';
W = F;
for iter = 1:length(s_in)
    X(K*J-K+1:K*J) = [];
    X = [s_in(iter, :).'; X];
    y(iter) = W.'*X;
    W = P*(W - mu*y(iter)*X) + F;
end


figure(666);
plot(s_in(:,1))
hold on
plot(y)
plot(t_sign_in(:,1))
ylim([-1 1]);
%% полученная ДН МР
phi = -90:0.1:89;
%steervec((n-1)*d/lambda, 10)  =  exp(1i*2*pi  /lambda*(n-1)*d  *sind(10) )

lambda = c/fsin;

b_n = zeros(1,length(phi));
b_c = zeros(1,length(phi));
Wn = ones(J, 1);
Wc = reshape(W,[Nel,J]).';
 for ii=1:length(phi)
     for n = 1:Nel
        s = steervec((n-1)*d/lambda,phi(ii));
        b_n(ii) = b_n(ii) + Wn'*s;                      % без коэффициентов
        b_c(ii) = b_c(ii) + Wc(:, n)'*s;                % с  коэффициентами
     end 
 end
F_n = abs(b_n).^2; 
F_c = abs(b_c).^2;

Fn_dB = 10*log10(F_n/max(F_n));
Fc_dB = 10*log10(F_c/max(F_c));

Fn_dB(Fn_dB<-40) = -40;
Fc_dB(Fc_dB<-40) = -40;

figure(5); 
grid on;
hold on;
% plot(phi,Fn_dB);
plot(phi,Fc_dB,'r'); 
yyy = ylim;
line([tetaSign tetaSign], yyy,'LineWidth',2)
line([tetaNois tetaNois], yyy,'LineWidth',2)







