clc;clear;close all;
load('NEWnoise_voice34mic_5.mat','VarName1');
fs = 8000;

X = VarName1;
voic1 = X(1:length(X)/3);
voic2 = X(length(X)/3+1:2*length(X)/3);
FROmk = X(2*length(X)/3+1:end);

NQWER = 160;
K_ys = 1024;

t_duration = length(voic1)/fs;
t = 0:1/fs:t_duration-1/fs;

Nel = 2;
c = 343;
ff = 2000;
d = c/ff/2;
% Frost....................................................................
J = 25;
K = Nel;
% mu = 0.0000005;
% mu = 2000000;
mu = 2^21;

sigArray(:,1) = voic1;
sigArray(:,2) = voic2;


W = zeros(K,J);
X = zeros(K,J);
y_true = zeros(length(sigArray),1);
if(J == 1)
    C = 1;
else
    C = (fir1(J-1,0.99,'low',chebwin(J,30)).'*K_ys).';
end

for ii = 1:length(W(:,1))
    W(ii,:) = (1/K*C.');
end

for iter = 1:length(sigArray(:,1))
    X(:,J) = [];
    X = [sigArray(iter, :).' X];
    y_true(iter) = sum(sum(W.*X));
    W = W - (y_true(iter)*X/mu');
    W = W + (1/K*(C-ones(1,K)*W));
end
y_true = (y_true/K_ys).';

W_mic1_lr = zeros(1,J);
W_mic2_lr = zeros(1,J);

X_mic1 = zeros(1,J);
X_mic2 = zeros(1,J);
y = zeros(1,length(sigArray));
if(J == 1)
    C = 1;
else
    C = floor(fir1(J-1,0.99,'low',chebwin(J,30)).'*K_ys).';
end  

C_lr = fliplr(C);
W_mic1_lr = floor(fliplr(1/K*C_lr));
W_mic2_lr = floor(fliplr(1/K*C_lr));
X_mic1 = sigArray(:,1).';
X_mic2 = sigArray(:,2).';
%%
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
% sum(y-y_try)

%%
% y_old = y;
% save('test_y.mat','y_old');

load('test_y.mat','y_old');
%%
figure(3)
plot(y);
hold on
plot(y_old);

% figure(4)
% plot(y_true);
% hold on
% plot(y);


% figure(5)
% plot(y_true);
% hold on
% plot(y-y_true);
% 
% figure(6)
% plot(FROmk/K_ys);
% hold on
% plot(y);
%%
hap = dsp.AudioPlayer('SampleRate',fs);
step(hap,voic1/400);
%%
hap = dsp.AudioPlayer('SampleRate',fs);
step(hap,y.'/400);
%%
step(hap,y_old.'/400);

