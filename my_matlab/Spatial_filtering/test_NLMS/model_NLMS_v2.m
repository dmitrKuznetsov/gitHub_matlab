% clc;clear;
close all;

load('mic1.txt');
mic1 = int32(reshape(mic1.',[],1));


load('ImResp_CompFIR.mat')
d = int32(floor(filter(bhi, 1, mic1)/(2^31)));

order = 32;
POSTSHIFT = 0;
Value_flip = flipud(Value);

mu = int32(2^31-1);

e = int32(zeros(length(d),1));
y = int32(zeros(length(d),1));
x = int32(zeros(order, 1));
b = int32(zeros(order, 1));

x0 = int32(0);
energy = int64(0);
acc = int64(0);
acc_l = int32(0);
acc_h = int32(0);
uShift = uint32(POSTSHIFT+1);
lShift = uint32(32) - uShift;
oneByEnergy = int32(0);
for k = 1:length(d)
    x = [mic1(k); x(1:end-1)];
    
    energy = bitshift  (bitshift(energy, 32) - bitshift(x0*x0, 1),-32);
    energy = bitshift  (bitshift(energy, 32) + bitshift(mic1(k)*mic1(k), 1),-32);
    acc = 0;
    for tapCnt = 1: order
        acc = acc + x(tapCnt)*b(tapCnt);
    end
    
    acc_l = bitand(acc, uint32(4294967295));
    acc_h = bitand(bitshift(acc, -32), uint32(4294967295));
    
    acc = bitor(bitshift(acc_l,-lShift),bitshift(acc_h,uShift));
    y(k) = int32(acc);
    e(k) = d(k) - y(k);
    
    
    
    
%     E = mod2n((x' * x ),32);
    E = (x' * x )/(2^31);
    b = b + floor(mu / E * x * e(k)'/(2^31));
end


figure();
hold on
%plot(bhi/max(bhi))
plot(b/max(b))
plot(Value/max(Value))



































