clc;clear;close all;
R = 4;          %���� �����
D = 2;          %��������
N = 4;          %������� �������
cicdec1 = dsp.CICDecimator(R,D,N); 
y_in = [1 2 3 4 5 zeros(1,16*5-5)];
% y_in = ones(1,100);
y_out = cicdec1(y_in.');


y_out_my = CIC_nmy_N4(y_in,R,38);
figure(101);
hold on
plot(y_out)
plot(y_out_my)










