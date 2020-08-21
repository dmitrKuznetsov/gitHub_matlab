function [] = plotSpectr(sig, N, numFig,fs)


Spec = abs(fft(sig,N));
f = (0:length(Spec)-1)/length(Spec)*fs;
figure(numFig)
plot(f,20*log10(Spec/max(Spec)) )
xlabel('Freq, Hz')
ylabel('Spectr')
end

