clc;clear;close all;
load('noise_voice34mic_4.mat','VarName1');
fs = 8000;
X = 80*VarName1/40000;
 
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
f = (-length(voic1)/2:length(voic1)/2-1)*fs/length(voic1);

s_in(:,1) = voic1;
s_in(:,2) = voic2;


IN_BUFF_SIZE = 64;
OUT_BUFF_SIZE = 8;
T_INTERV = 20;
NUM_INTERV = 20/T_INTERV;


Nel = 2;
c = 343;
ff = 2000;
d = c/ff/2;
%% ...................................Frost..................... ............
J = 32;
K = Nel;
K_ys = 1024;
NQWER = 160;

mu = 2^21;

W = zeros(K,J);
X = zeros(K,J);
y_true = zeros(length(s_in),1);
if(J == 1)
    C = 1;
else
    C = (fir1(J-1,0.99,'low',chebwin(J,30)).'*K_ys).';
end

for ii = 1:length(W(:,1))
    W(ii,:) = (1/K*C.');
end

for iter = 1:length(s_in(:,1))
    X(:,J) = [];
    X = [s_in(iter, :).' X];
    y_true(iter) = sum(sum(W.*X));
    W = W - (y_true(iter)*X/mu');
    W = W + (1/K*(C-ones(1,K)*W));
end
y_true = (y_true/K_ys).';
%% ................................... MK ...................................
W_mic1_lr = zeros(1,J);
W_mic2_lr = zeros(1,J);

X_mic1 = zeros(1,J);
X_mic2 = zeros(1,J);
y = zeros(1,length(s_in));
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
% sum(y-y_try)



figure(4)
hold on
plot(y_true/max(y_true));
plot(y/max(y));




