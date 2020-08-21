clc;clear;
load('noise_voice34mic_4.mat','VarName1');
close all;
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
mu = 2000000;
sigArray(:,1) = voic1;
sigArray(:,2) = voic2;

% X = zeros(K*J,1);
% y = zeros(length(sigArray(:,1)),1);
% C = zeros(K*J,J);
% W = zeros(K*J,1);
% I = diag(ones(1,K*J));
% P = zeros(K*J,K*J);
% F = zeros(K*J,1);                                    
% for j = 1:J
%     C(:,j) = [zeros(1,(j-1)*K) ones(1,K) zeros(1,J*K-j*K)].';
% end
% if(J == 1)
%     FF = 1;
% else
%     FF = round(fir1(J-1,0.99,'low',chebwin(J,30)).'*10^4);
% end
% % figure(111);
% % freqz(FF,1)
% F = C/(C.'*C)*FF;                           %!!!C* inv(C.'*C)*FF;
% P = I - C/(C.'*C)*C.';
% W = F;
% for iter = 1:length(sigArray)
%     X(K*J-K+1:K*J) = [];
%     X = [sigArray(iter, :).'; X];
%     y(iter) = W.'*X;
%     W = P*(W - mu*y(iter)*X/mu) + F;
%     W = round(W);
% end
%%
W = zeros(K,J);
X = zeros(K,J);
y = zeros(length(sigArray),1);
if(J == 1)
    C = 1;
else
    C = floor(fir1(J-1,0.99,'low',chebwin(J,30)).'*10^3);
end

for ii = 1:length(W(:,1))
    W(ii,:) = 1/K*C.';
end

for iter = 1:length(sigArray)
    X(:,J) = [];
    X = [sigArray(iter, :).' X];
    y(iter) = sum(sum(W.*X));
    W = W - y(iter)*X/mu;
    W = W + 1/K*(C.'-ones(1,K)*W);
    W = floor(W);
end



y = floor(y/10^3);
figure(3)
plot(t,voic2);
hold on
plot(t,y);
xlabel('Time (sec)'); ylabel ('Amplitude');
title('Frost Beamformer Output'); 








% %%
% hap = dsp.AudioPlayer('SampleRate',fs);
% step(hap,voic1/400);
%  %%
% hap = dsp.AudioPlayer('SampleRate',fs);
% step(hap,y/400);








