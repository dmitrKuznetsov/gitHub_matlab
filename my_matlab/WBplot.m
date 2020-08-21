function WBplot(W, doa, Ts, varargin)
% WBPLOT
% Required Input Arguments:
% - W : J x M weight matrix, where: M = number of antennas
% J = number of taps or length
% of the temporal filter
% - doa : direction of arrival of each mobile user (SOI's DoA
% must be the first element in the array);
% - Ts : Sample time
%
% Additiona Input Argument
% - Nf : number of frequency bins. Default = 50
% - Na : number of angle bins. Default = 60
% - Dim : Graph dimension ['2D','3D']. Default = 3D
% - Color : Colormap ['jet','white','gray']. Default = jet
% - d : Normalized spacing between sensors. Default = 0.5
% - c : Wavespeed. Default = 3E8
%%%%%%%%%%%%%%%%%%%% Input Parameter Validation %%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
defaultDim = '3D';
defaultNf = 50;
defaultNa = 60;
defaultd = 0.5;
defaultc = 3E8;
defaultColor = 'jet';
expectedDim = {'2D','3D'};
expectedColor = {'white','gray','jet'};
addRequired(p,'w',@isnumeric);
addRequired(p,'doa',@isnumeric);
addRequired(p,'Ts',@isnumeric);
addOptional(p,'Dim',defaultDim,@(x) any(validatestring(x,expectedDim)));
addOptional(p,'Nf',defaultNf,@isnumeric);
addOptional(p,'Na',defaultNa,@isnumeric);
addOptional(p,'d',defaultd,@isnumeric);
addOptional(p,'c',defaultc,@isnumeric);
addOptional(p,'Color',defaultColor,@(x) any(validatestring(x,expectedColor)));
parse(p, W, doa, Ts, varargin{:});
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = 1;
K = length(doa);
[J M] = size(W);
angRes = p.Results.Na - 1;
freqRes = p.Results.Nf - 1;
response = zeros(angRes + 1, freqRes + 1);
fs = 1/p.Results.Ts;
mu = p.Results.d/(p.Results.c * p.Results.Ts);
mu = 1;
omega = 2*pi*p.Results.Ts; % This should be multiplied by freq
ian = 1;
ifr = 1;
temp = 0;
for ang = -pi/2:pi/angRes:pi/2
ifr = 1;
for freq = 0:(fs/2)/freqRes:fs/2
temp = 0;
for countj = 0:J-1
temp = temp + (exp(-1j*omega*freq*((0:M-1)*mu*sin(ang) + countj)))*(W(countj+1,:)');
end
response(ian,ifr) = temp;
ifr=ifr+1;
end
ian=ian+1;
end
ang=-90:180/angRes:90;
freq = 0:1/freqRes:1;
if strcmp(p.Results.Dim,'3D')
figure
[x y] = meshgrid(freq,ang);
surfc(x, y,10*log10(abs(response)));
colormap(p.Results.Color);
% surfc(x, y,10*log10(abs(response)),'EdgeColor','none'); % This can be used so the grids on the surface are not w
% zlim([-40 10]);
ylim([-90 90]);
xlabel('Normalized Frequency (\Omega/\pi)');
ylabel('DOA (\theta)');
zlabel('Gain (dB)');
view(60, 20);
end
if strcmp(p.Results.Dim,'2D')
figure
for k = 1:1:freqRes+1
temp = 10*log10(abs(response(:,k)));
plot(ang,temp,'k')
hold on
end
xlim([-90 90]);
xlabel('DOA (\theta)');
ylabel('Gain (dB)');
end
end