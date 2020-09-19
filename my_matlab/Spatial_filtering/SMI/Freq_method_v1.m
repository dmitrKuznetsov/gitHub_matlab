clc; clear; close all;
c = 343;
Nel = 2;

f_required = 2000;
d = c/f_required/2;

fs = 8e3;
t_duration = 2;
t = 0:1/fs:t_duration-1/fs;
NSampPerFrame = length(t);

%% .................................Исходные Сигналы...............................
[signl,~] = audioread('dft_voice_8kHz.wav');
[noise,~] = audioread('nois.wav');
if(length(signl)>=NSampPerFrame || length(noise)>=NSampPerFrame)
    signl(NSampPerFrame+1:end) = [];
    noise(NSampPerFrame+1:end) = [];
else
    warning('проблема с размером сигналов');
end

FFsign = fft(signl);
FFsign(1) = 0;
signl = ifft(FFsign);

FFnois = fft(noise);
FFnois(1) = 0;
noise = ifft(FFnois);

signPwr = mean((signl).^2);
noisPwr = mean((noise).^2);

SNR_true = signPwr/noisPwr;
SNR_true_db = pow2db(SNR_true);
%% ..................................Принятые сигналы...............................
tetaSign = 0;                    
tetaNois = 55;

t_sign_in = Received_from_direction(signl, Nel, tetaSign, d, fs);
t_nois_in = Received_from_direction(noise, Nel, tetaNois, d, fs);
s_in = t_sign_in + t_nois_in;

% for iii = 1: Nel
%     s_in(:,iii) = s_in(:,iii) + sqrt(noisPwr/1000)*randn(length(s_in),1);
% end

%% ....................................Алгоритм...................................

X1 = fft(s_in(:, 1));
X2 = fft(s_in(:, 2));

X1_pos = X1(1:length(X1)/2);
X2_pos = X2(1:length(X2)/2);

L = length(X1_pos)/64;
K = length(X1_pos)/L;


% imp = fir1(2*L-1,0.99,'low',chebwin(2*L,30)).';
imp = [1; zeros(L-1,1)];
freq_resp = fft(  imp, length(imp)*2  );
F = freq_resp(1:length(freq_resp)/2);

C = ones(2,1);




index = zeros(L, K);
X_array = zeros(2, K);
W_array = zeros(2, L);
for i = 0:L-1
    index(i+1,:) = i*K+1 : (i+1)*K;
    X_array(1,:) = X1_pos(index(i+1,:));
    X_array(2,:) = X2_pos(index(i+1,:));
    R_ll = 1/K*(X_array*X_array.');
    R_ll_inv = inv(R_ll);

    W_opt(:,i+1) = R_ll_inv*C*  inv(C.'*R_ll_inv*C)  *F(i+1);
end 




%% ....................................Построение ДН от частоты...................................

  
freqbin = linspace(0,fs/2-fs/2/L, L);                                   
phi = -90:0.5:90-0.5;

BP = zeros(length(freqbin), length(phi));

for kk = 1:length(freqbin)
    %steervec((n-1)*d/lambda, 10)'  =  exp(-1i*2*pi  /lambda*(n-1)*d  *sind(10) )
    lambda = c/freqbin(kk);
    P = zeros(1,length(phi));
    for ii=1:length(phi)
        stv = steervec(( (1:Nel) -1)*d/lambda, phi(ii))';
        P(ii) =  stv*W_opt(:,kk);                
    end
    BP(kk,:) = abs(P).^2;
end
BP_dB = 10*log10(BP);

BP_dB(BP_dB<-40) = -40;
%%
figure(5); 
grid on;
hold on;
plot(phi,BP_dB,'r'); 
yyy = ylim;
line([tetaSign tetaSign], yyy,'LineWidth',1,'Color','k')
line([tetaNois tetaNois], yyy,'LineWidth',1)



 
[PHI,FRREQ] = meshgrid(phi,freqbin);

figure(777);
colormap(gray);
hold on;
grid on;
surf(PHI,FRREQ,BP_dB,'LineStyle','none'); 

xlabel('Угол прихода, град.');
ylabel('Частота, Гц');
ylim([0 3800])





