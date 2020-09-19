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

X1 = fft(s_in(:,1));
X2 = fft(s_in(:,2));


L = 16*100;
Rxx = zeros(2*L,2*L);

f = (0:length(s_in(:,1))-1)/t_duration;
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







%% ....................................Построение ДН от частоты...................................
num_bin = 2048;
wfreq = fft(reshape(W,Nel,J).',num_bin);                                    % из импульсной характеристики каналов получаем частотную
freqvec = [0:num_bin/2 -num_bin/2+1:-1]*fs/num_bin;                         % частоты, которые соответствуют отсчетам частотной характеристики wfreq

freqbin = freqvec(1  :  length(freqvec)/2 - length(freqvec)/2/10);            % частоты, для которых строю ДН
phi = -90:0.5:90-0.5;

BP = zeros(length(freqbin), length(phi));

for kk = 1:length(freqbin)
%     bandidx = find( freqvec == freqbin(kk)) ;                                 % индекс freqvec, который соответ частоте freqbin(kk)

    %steervec((n-1)*d/lambda, 10)'  =  exp(-1i*2*pi  /lambda*(n-1)*d  *sind(10) )
    lambda = c/freqbin(kk);

    P = zeros( 1, length(phi) );

    for ii = 1:length(phi)
        stv = steervec(( (1:Nel) -1)*d/lambda,phi(ii))';                  
        P(ii) =  stv * (   wfreq(kk, :).'   );                
    end
    BP(kk,:) = abs(P).^2;
end

BP_dB = 10*log10(BP);

BP_dB(BP_dB < -40) = -40;



num_bin = 2*L;
freqvec = [0:num_bin/2 -num_bin/2+1:-1]*fs/num_bin;         % частоты, которые соответствуют отсчетам частотной характеристики wfreq
NN = L;
freqbin = freqvec(1:length(freqvec)/2/NN:length(freqvec)/2-length(freqvec)/2/10);                                     % частоты, для которых строю ДН
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
        b_n(ii) =  stv*conj(w_rm);                      
        b_c(ii) =  stv*conj(W_opt(:,bandidx));                
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



