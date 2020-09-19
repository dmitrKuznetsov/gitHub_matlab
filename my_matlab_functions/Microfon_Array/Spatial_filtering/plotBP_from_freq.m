function [] = plotBP_from_freq(W, Nel, J, tetaSign, tetaNois)
c = 343;
f_required = 2000;
lambda = c/f_required;
d = lambda/2;
fs = 8e3;

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
%%
figure();
grid on;
hold on;
plot(phi,BP_dB,'LineWidth',0.25,'Color',[0.501960813999176 0.501960813999176 0.501960813999176]); 
yyy = ylim;
line([tetaSign tetaSign], yyy,'LineWidth',1,'Color','k')
line([tetaNois tetaNois], yyy,'LineWidth',1,'Color','k')
ylabel('BP, дБ')
xlabel('Угол прихода, град')
xlim([-90 90])


 %%
[PHI,FRREQ] = meshgrid(phi, freqbin);

figure1 = figure();
colormap(gray);

axes1 = axes('Parent',figure1);
hold(axes1,'on');


surf(PHI,FRREQ,BP_dB,'Parent',axes1,'LineStyle','none'); 

ylabel('Частота, Гц','Rotation',-80)
xlabel('Угол прихода, град')
zlabel('BP, дБ')
xlim([-90 90])
ylim([1 freqbin(end)])
view(axes1,[174.34 62]);
grid(axes1,'on');

% set(axes1,'FontSize',14);
% zlim([-40 5])

end

