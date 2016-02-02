%% Parameters
fc = 2e9; % carry frequency (Hz)
ts = 1e-4; % sample interval (s)
v = 25/3; % transmitter-receiver relative speed (m/s)
c = 3e8; % speed of light (m/s)
fd = v*fc/c; % maximum Doppler shift
Ns = 10000; % number of samples in simulation
x = (randn(Ns, 1) + sqrt(-1)*randn(Ns, 1))/sqrt(2); % Gaussian noise

%% The filter method
N = 1000; % the length of window 2*N+1
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
c = conv(x,g_hat); % channel gain
c = c(2*N+1:end); % discard samples due to filter transients

%% Simulation task
%-----fading envelope-----%
figure
plot(abs(c));
xlabel('time samples')
ylabel('channel gain')
title('Channel Gain')

%-----pdf-----%
figure
xx = 0:0.01:4;
pf = raylpdf(xx,1/sqrt(2));
plot(xx,pf); grid on; hold on
xlabel('x')
ylabel('P(x)')
ksdensity(abs(c),'support','positive');
legend('theory','empirical')
title('Empirical PDF')


%-----cdf-----%
figure
cdfplot(abs(c));hold on
cf = raylcdf(xx,1/sqrt(2));
plot(xx,cf);
legend('theory','empirical')

%-----autocorrelation-----%
figure
autocorr(abs(c));
tt = 0:ts:20*ts;
J_autocorr = besselj(0,2*pi*fd*tt);
hold on
plot(0:20,J_autocorr);
legend('empirical','theory')

%-----psd-----%
figure
periodogram(c);
f = -fd:fd/Ns:fd;
Sc = 1./(pi*fd*sqrt(1-(f/fd).^2)); % Doppler spectrum
hold on
ff = f/(floor(fd))+1;
plot(ff,10*log10(Sc),'r');
legend('empirical','theory')
