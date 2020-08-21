clc;clear;close all;
c = 343;
Nel = 4;

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
NoisdftFileReader = dsp.AudioFileReader('FemaleSpeech-16-8-mono-3secs.wav',...%FemaleSpeech-16-8-mono-3secs.wav   nois.wav   my_nois  'Laughter-16-8-mono-4secs.wav'
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
[y_LMS,~] = my_Frost(s_in, J, K, mu);

[y_RLS, W] = func_LC_RLS_broadband(s_in,K,J,1);
%%
figure(1)
hold on
% plot(s_in(:,1))
plot(y_LMS)
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
num_bin = 32;%length(sig);
wfreq = fft(reshape(W,Nel,[]).',num_bin);
freqvec = [0:num_bin/2 -num_bin/2+1:-1]*fs/num_bin;
freqbin = 250:250:4000;
for kk = 1:length(freqbin)
    kk;
    bandidx = find( freqvec==freqbin(kk)) ;  % corresponding to 1kHz
    freqvec(bandidx);

    % WBplot(reshape(W,Nele,[]).', 0, 1/fs,'Nf',128)

    phi = -90:0.5:90-0.5;
    %steervec((n-1)*d/lambda, 10)  =  exp(1i*2*pi  /lambda*(n-1)*d  *sind(10) )

    lambda = c/freqbin(kk);

    b_n = zeros(1,length(phi));
    b_c = zeros(1,length(phi));

    w_rm = ones(Nel,1);
    for ii=1:length(phi)
        stv = steervec(( (1:Nel) -1)*d/lambda,phi(ii)).';
        b_n(ii) =  stv*w_rm;                      
        b_c(ii) =  stv*(wfreq(bandidx,:)).';                
    end
    F_n(kk,:) = abs(b_n).^2; 
    F_c(kk,:) = abs(b_c).^2;
end
Fn_dB = 10*log10(F_n);
Fc_dB = 10*log10(F_c);

Fn_dB(Fn_dB<-40) = -40;
Fc_dB(Fc_dB<-40) = -40;
%%
figure(5); 
grid on;
hold on;
% plot(phi,Fn_dB);
plot(phi,Fc_dB,'r'); 
yyy = ylim;
line([tetaSign tetaSign], yyy,'LineWidth',1,'Color','k')
line([tetaNois tetaNois], yyy,'LineWidth',1)

% legend('Fn dB','Fc dB');

 
[PHI,FRREQ] = meshgrid(phi,freqbin);


figure1 = figure(777);
colormap(gray);
% Create axes
axes1 = axes('Parent',figure1);
hold(axes1,'on');

% Create surf
surf(PHI,FRREQ,Fc_dB,'Parent',axes1); 
% Create ylabel
ylabel('Частота, Гц');

% Create xlabel
xlabel('Угол прихода, град.');
ylim([0 3800])
view(axes1,[177.22 31.92]);
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16);
% zlim([-40 5])






