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
        F_sign(index,i) = FFsign(index)*exp( -1i*2*pi*f_00/c*d*(i-1)*sind(tetaSign) );
        F_nois(index,i) = FFnois(index)*exp( -1i*2*pi*f_00/c*d*(i-1)*sind(tetaNois) );
    end
end

F_sign(end) = real(F_sign(end));
F_nois(end) = real(F_nois(end));

for i = 1:Nel
    F_sign(pos_len+1:2*pos_len-2,i) = conj(F_sign(pos_len-1:-1:2,i));    
    F_nois(pos_len+1:2*pos_len-2,i) = conj(F_nois(pos_len-1:-1:2,i));
end



t_sign_in = ifft(F_sign,[]);
t_nois_in = ifft(F_nois,[]);
s_in = t_sign_in + t_nois_in;
%%
X1 = fft(s_in(:,1));
X2 = fft(s_in(:,2));
L = 16*100;
Rxx = zeros(2*L,2*L);

f = (0:length(s_in(:,1))-1)/t_duration;
% Обрабатываю положительную часть спектра
f_pos = f(2:length(f)/2+1);                                       

                                                
freq_per_bin = max(f_pos)/L;                                
sample_per_bin = length(f_pos)/L;                           
%......................................................................
C = ones(2,1);

% F = 20*ones(L,1);
imp = fir1(2*L-1,0.99,'low',chebwin(2*L,30)).';
freq_resp = fft(imp);
F = freq_resp(2:length(freq_resp)/2+1);
%......................................................................

index = zeros(L, sample_per_bin);
X_array = zeros(2, sample_per_bin);
W_array = zeros(2, L);
for i = 0:L-1
    index(i+1,:) = i*sample_per_bin+2:(i+1)*sample_per_bin+1;
    X_array(1,:) = X1(index(i+1,:));
    X_array(2,:) = X2(index(i+1,:));
    R_ll = 1/sample_per_bin*(X_array*X_array.');
    R_ll_inv = inv(R_ll);

    W_opt(:,i+1) = R_ll_inv*C*inv(C.'*R_ll_inv*C)*F(i+1);
end 

%%
num_bin = 2*L;
freqvec = [0:num_bin/2 -num_bin/2+1:-1]*fs/num_bin;         % частоты, которые соответствуют отсчетам частотной характеристики wfreq
NN = L;
freqbin = freqvec(1 : length(freqvec)/2/NN : length(freqvec)/2-length(freqvec)/2/10);                                     % частоты, для которых строю ДН
% freqbin = 250:250:4000;                                     % частоты, для которых строю ДН
phi = -90:0.5:90-0.5;
F_n = zeros(length(freqbin),length(phi));
F_c = zeros(length(freqbin),length(phi));

for kk = 1:length(freqbin)
    bandidx = find( freqvec==freqbin(kk)) ;                 % индекс freqvec, который соответ частоте freqbin(kk)
    %steervec((n-1)*d/lambda, 10)'  =  exp(-1i*2*pi  /lambda*(n-1)*d  *sind(10) )
    lambda = c/freqbin(kk);

    b_n = zeros(1,length(phi));
    b_c = zeros(1,length(phi));

    w_rm = ones(Nel,1);
    for ii=1:length(phi)
        stv = steervec(( (1:Nel) -1)*d/lambda,phi(ii))';
        b_n(ii) =  stv*(w_rm);                      
        b_c(ii) =  stv*(W_opt(:,bandidx));                
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
surf(PHI,FRREQ,Fc_dB,'Parent',axes1,'LineStyle','none'); 
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







































