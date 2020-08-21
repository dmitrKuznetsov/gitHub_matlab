clc;clear;close all;
c = 343;
Nel = 8;

f_required = 2000;
lambda = c/f_required;
d = lambda/2;

fs = 8e3;
t_duration = 1.5;
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);
 
% hap = dsp.AudioPlayer('SampleRate',fs);
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

for i = 1:Nel
F_sign(pos_len+1:2*pos_len-2,i) = conj(F_sign(pos_len-1:-1:2,i));    
F_nois(pos_len+1:2*pos_len-2,i) = conj(F_nois(pos_len-1:-1:2,i));
end

t_sign_in = ifft(F_sign,[],1,'symmetric');
t_nois_in = ifft(F_nois,[],1,'symmetric');
s_in = t_sign_in + t_nois_in;
% ................................Обработка.............................
J = 30;
K = Nel;
mu = 0.005;
% [y_LMS,~] = my_Frost(s_in, J, K, mu);

[y_RLS, W] = func_LC_RLS_broadband(s_in,K,J,1);
%%
figure(1)
hold on
% plot(s_in(:,1))
% plot(y_LMS)
plot(y_RLS)
plot(t_sign_in(:,1))
legend('LMS',"RLS",' речь')
% figure(666);
% plot(s_in(:,1))
% hold on
% plot(y)
% plot(t_sign_in(:,1))
% ylim([-1 1]);


%%
num_bin = 1000;
wfreq = fft(reshape(W,Nel,[]).',num_bin);
freqvec_pos = (0:num_bin/2)*fs/num_bin;

% phi = tetaSign;
phi = tetaNois;




%steervec((n-1)*d/lambda, 10)  =  exp(1i*2*pi  /lambda*(n-1)*d  *sind(10) )
b_n = zeros(1,length(freqvec_pos));
b_c = zeros(1,length(freqvec_pos));
w_rm = ones(Nel,1);
for kk = 1:length(freqvec_pos)
    lambda = c/freqvec_pos(kk);
    stv = steervec(( (1:Nel) -1)*d/lambda,phi).';
    b_n(kk) =  stv*w_rm;                      
    b_c(kk) =  stv*(wfreq(kk,:)).';                
end
F_n = abs(b_n).^2; 
F_c = abs(b_c).^2;

Fn_dB = 10*log10(F_n/max(F_n));
% Fc_dB = 10*log10(F_c/max(F_c));
Fc_dB = 10*log10(F_c);


figure1 = figure;

% Create axes
axes1 = axes('Parent',figure1); 
hold(axes1,'on');

% Create multiple lines using matrix input to plot

plot(freqvec_pos,Fn_dB,'Parent',axes1,'LineWidth',4,'Color',[0 0 0]);
plot(freqvec_pos,Fc_dB,'Parent',axes1,'LineWidth',4,'Color',[0 0 0],'LineStyle','-.'); 

% % Create ylabel
% ylabel('BP, дБ');
% 
% % Create xlabel
% xlabel('Частота, Гц');
ylim([-60 0])

grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16);







