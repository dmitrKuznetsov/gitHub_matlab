clc;clear;close all;
load('noise_voice34mic_4.mat','VarName1');
fs = 8000;

X = VarName1;
voic1 = X(1:length(X)/2);
voic2 = X(length(X)/2+1:end);
% voic1(1:400) = [];
% voic2(1:400) = [];
for kkk = 1:length(voic1)
    if(voic1(kkk)>1024)
        voic1(kkk)=1024;
    end   
    if(voic2(kkk)>1024)
        voic2(kkk)=1024;
    end 
    if(voic1(kkk)<-1024)
        voic1(kkk)=-1024;
    end 
    if(voic2(kkk)<-1024)
        voic2(kkk)=-1024;
    end 
end

NQWER = 160;
K_ys = 1024;



% voic1(NQWER+1:end) = [];
% voic2(NQWER+1:end) = [];



t_duration = length(voic1)/fs;
t = 0:1/fs:t_duration-1/fs;
% figure(1);
% plot(t, voic1);
% hold on
% plot(t, voic2);
% ylim([-512 512])


Nel = 2;
c = 343;
ff = 2000;
d = c/ff/2;
% Frost....................................................................
J = 5;
K = Nel;
% mu = 0.0000005;
% mu = 2000000;
mu = 2^21;

sigArray(:,1) = voic1;
sigArray(:,2) = voic2;

W = zeros(K,J);
X = zeros(K,J);
y_try = zeros(length(sigArray),1);
if(J == 1)
    C = 1;
else
    C = fix(fir1(J-1,0.99,'low',chebwin(J,30)).'*K_ys).';
end

for ii = 1:length(W(:,1))
    W(ii,:) = fix(1/K*C.');
end

for iter = 1:length(sigArray(:,1))
    X(:,J) = [];
    X = [sigArray(iter, :).' X];
    y_try(iter) = sum(sum(W.*X));
    W = W - fix(y_try(iter)*X/mu');
    W = W + fix(1/K*(C-ones(1,K)*W));
%     W = fix(W)

end
y_try = fix(y_try/K_ys).';

%%


W_mic1_lr = zeros(1,J);
W_mic2_lr = zeros(1,J);

X_mic1 = zeros(1,J);
X_mic2 = zeros(1,J);
y = zeros(1,length(sigArray));
if(J == 1)
    C = 1;
else
    C = fix(fir1(J-1,0.99,'low',chebwin(J,30)).'*K_ys).';
end  

C_lr = fix(C);
W_mic1_lr = (fliplr(1/K*C_lr));
W_mic2_lr = (fliplr(1/K*C_lr));
X_mic1 = sigArray(:,1).';
X_mic2 = sigArray(:,2).';

Niter = floor(length(X_mic1)/160);
data_mic1 = zeros(1, (J-1)+160);
data_mic2 = zeros(1, (J-1)+160);
%%


for i_iter = 0:Niter-1
    index = 1+NQWER*i_iter : NQWER+NQWER*i_iter;
    data_mic1(J:J+NQWER-1) = X_mic1(index);
    data_mic2(J:J+NQWER-1) = X_mic2(index);
    for i_elem = 1:NQWER %!!!!!!!
        y(i_elem+NQWER*i_iter) = 0;
        for i_otvod = 0:J-1
            y(i_elem+NQWER*i_iter) = y(i_elem+NQWER*i_iter) + data_mic1(i_elem+i_otvod)*W_mic1_lr(i_otvod+1)...
                                                            + data_mic2(i_elem+i_otvod)*W_mic2_lr(i_otvod+1);
            y(i_elem+NQWER*i_iter) = y(i_elem+NQWER*i_iter);
        end
        for i_otvod = 1:J%!!!
            W_mic1_lr(i_otvod) = W_mic1_lr(i_otvod) - y(i_elem+NQWER*i_iter)*data_mic1(i_elem+(i_otvod-1))/mu;
            W_mic2_lr(i_otvod) = W_mic2_lr(i_otvod) - y(i_elem+NQWER*i_iter)*data_mic2(i_elem+(i_otvod-1))/mu;
        end
        for i_otvod = 1:J%!!!
            V_corect(i_otvod) = (C_lr(i_otvod) - W_mic1_lr(i_otvod) - W_mic2_lr(i_otvod))/K;
            W_mic1_lr(i_otvod) = W_mic1_lr(i_otvod) + V_corect(i_otvod);
            W_mic2_lr(i_otvod) = W_mic2_lr(i_otvod) + V_corect(i_otvod);
        end
    

    end
    
    %в преамбулу добавим данные из прошлой итерации
    data_mic1(1:J-1) = data_mic1(end-J+2:end);
    data_mic2(1:J-1) = data_mic2(end-J+2:end);
        
end        
y = y/K_ys;
sum(y-y_try)


figure(3)
plot(t,voic2);
hold on
plot(t,y);
plot(t,y_try);
xlabel('Time (sec)'); ylabel ('Amplitude');
title('Frost Beamformer Output');
ylim([-400 400])

figure(4)
plot(t,voic2);
hold on
plot(t,y-y_try);
xlabel('Time (sec)'); ylabel ('Amplitude');
title('Ошибка из-за округления');
ylim([-400 400])
% ошибка обратно пропорционально K_ys - т.е увеличиваем число разрядов, 
% учавствующих в подсчете, если сигнал мал - ошибка влияет сильно







% %%
% hap = dsp.AudioPlayer('SampleRate',fs);
% step(hap,voic2/500);
% 
% %%
% hap = dsp.AudioPlayer('SampleRate',fs);
% step(hap,y.'/500);
%  %%
% hap = dsp.AudioPlayer('SampleRate',fs);
% step(hap,y_try.'/500);








