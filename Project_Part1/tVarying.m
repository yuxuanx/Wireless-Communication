%% Parameters
fc = 2e9; % carry frequency (Hz)
ts = 1e-4; % sample interval (s)
v = 25/3; % transmitter-receiver relative speed (m/s)
c = 3e8; % speed of light (m/s)
fd = v*fc/c; % maximum Doppler shift
Ns = 300; % number of samples in simulation
x = (randn(Ns, 1) + sqrt(-1)*randn(Ns, 1))/sqrt(2); % Gaussian noise

%% The filter method
N = 5; % the length of window 2*N+1
% window shape/@rectwin @hamming etc.
w = window(@rectwin,2*N+1);

%-----impulse response-----%
t = -N*ts:ts:N*ts; % make impulse response causal
Z = 2*pi*fd*abs(t)'; % argument of Bessel function
J = besselj(1/4,Z); % Bessel function of the first kind, order 1/4x
g = J./(abs(t)'.^(1/4));
g(N+1) = (pi*fd)^(1/4)/gamma(5/4);
g_hat = g.*w./sqrt(sum(abs(g.*w).^2)); % normalize impulse response

%-----channel gain-----%
c1 = conv(x,g_hat); % channel gain
c1 = c1(2*N+1:end); % discard samples due to filter transients

%% The spectrum method
fs = 1/ts; % sampling frequency
f_limit = floor(fd*Ns/fs); % |f| <= fd
f = 0:fs/Ns:f_limit*fs/Ns;
Sc = 1./(pi*fd*sqrt(1-(f/fd).^2)); % Doppler spectrum

%-----impulse response-----%
G = sqrt(Sc);

% for G_hat, when 0<f<fd G_hat(f) = G(f), when fs-fd<f<fs G_hat = G(f-fs)
G_hat = zeros(Ns,1);
G_hat(1:f_limit+1) = G(1:end);
G_hat(end-f_limit+1:end) = G(end:-1:2);

%-----channel gain-----%
X = x*sqrt(Ns)/std(G_hat); % make c2 have unit variance
c2 = ifft(G_hat.*X);

%% Simulation task
%-----fading envelope-----%
figure
plot(abs(c1));hold on
plot(abs(c2));

%-----pdf-----%
figure
ksdensity(abs(c1),'support','positive');hold on
ksdensity(abs(c2),'support','positive');hold on
xx = 0:0.01:4;
pf = raylpdf(xx,1/sqrt(2));
plot(xx,pf);

%-----cdf-----%
figure
cdfplot(abs(c1));hold on
cdfplot(abs(c2));hold on
cf = raylcdf(xx,1/sqrt(2));
plot(xx,cf);

%-----autocorrelation-----%
figure
subplot(2,1,1);
autocorr(abs(c1));
subplot(2,1,2);
autocorr(abs(c2));
tt = 0:ts:20*ts;
J_autocorr = besselj(0,2*pi*fd*tt);
figure
plot(0:20,J_autocorr);

%-----psd-----%
figure
subplot(2,1,1);
periodogram(c1);
subplot(2,1,2);
periodogram(c2);
f = -fd:fd/Ns:fd;
Sc = 1./(pi*fd*sqrt(1-(f/fd).^2)); % Doppler spectrum
figure
plot(f,10*log10(Sc));

