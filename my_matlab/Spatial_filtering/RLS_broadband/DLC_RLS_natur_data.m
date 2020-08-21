clc;clear;close all;
load('noise_voice34mic_4.mat','VarName1');
fs = 8000;
X = 80*VarName1/40000;
% X(1:300) = []; 
% voic1 = X; 
voic1 = X(1:length(X)/2);
voic2 = X(length(X)/2+1:end);
voic1(1:4000) = [];
voic2(1:4000) = [];

FF1 = fft(voic1);
FF2 = fft(voic2);
FF1(1) = 0;
FF2(1) = 0;
voic1 = ifft(FF1);
voic2 = ifft(FF2);

t_duration = length(voic1)/fs;
t = 0:1/fs:t_duration-1/fs;
% figure(1);
% plot(t, voic1);
% hold on
% plot(t, voic2);
% ylim([-0.8 0.8])


f = (-length(voic1)/2:length(voic1)/2-1)*fs/length(voic1);

Nel = 2;
c = 343;
ff = 2000;
d = c/ff/2;
hap = dsp.AudioPlayer('SampleRate',fs);



J = 35;
K = Nel;


sigArray(:,1) = voic1;
sigArray(:,2) = voic2;
mu = 0.1;
[my_FrostOut,W] = my_Frost(sigArray, J, K, mu);
[y_RLS, ~] = func_LC_RLS_broadband(sigArray,K,J,1);

figure(3)
hold on
plot(t,voic2);
plot(t,my_FrostOut);
plot(t,y_RLS);
xlabel('Time (sec)'); ylabel ('Amplitude (V)');


%% с одного из микрофонов voic1 - ближе к микро, voic2 - дальше от микро
step(hap,voic2);
%% мой LMS
release(hap);
step(hap,my_FrostOut);
%% мой RLS
release(hap);
step(hap,y_RLS);
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
% line([tetaSign tetaSign], yyy,'LineWidth',1,'Color','k')
% line([tetaNois tetaNois], yyy,'LineWidth',1)

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
ylabel('„астота, √ц');

% Create xlabel
xlabel('”гол прихода, град.');
ylim([0 3800])
view(axes1,[177.22 31.92]);
grid(axes1,'on');
% Set the remaining axes properties
set(axes1,'FontSize',16);
% zlim([-40 5])





















