clear all
close all
clc

% parameters
v  = 30*1000/3600; % relative speed
fc = 2e9;          % carrier frequency
c  = 3e8;
Ts = 10^(-4);      % Sample interval
fs = 1/Ts;
N  = 1000;         % No. of filter taps
Ns = 300;         % No. of samples

% Compute Doppler Frequency
fd = (v*fc)/c;

% Generate g(t) using Eqn.8
g = [besselj(0.25, 2*pi*fd*Ts*abs(-N:-1)')./(abs(-N:-1)'*Ts).^(1/4);
       (fd*pi)^(1/4)/gamma(5/4); 
       besselj(0.25, 2*pi*fd*Ts*(1:N)')./((1:N)'*Ts).^(1/4)];
g = g/norm(g);

G = fft(g,Ns); 

% Generate X(k)
X_k = sqrt(Ns)*(randn(1,Ns) +1i*randn(1,Ns));

% Compute c(t)
C_k = X_k.*G';
c_t = ifft(C_k);

% Plot: Fading envelope
figure(1);
plot(abs(c_t));

title('Fading envelope of c(nTs)');
xlabel('t');
ylabel('Amplitude');
legend('Spectrum Method');

%Plot: CDF
h_t = hist(abs(c_t),100);
cdf = cumsum(h_t) / Ns;
figure(2)
plot(cdf);
legend('Spectrum Method');
title('Cumulative Distribution Function');
xlabel('c');
ylabel('cdf(c)');

% Plot: PDF
pdf = diff(cdf);
figure(3)
plot(pdf);
legend('Spectrum Method');
title('Probability Density Function');
xlabel('c');
ylabel('pdf(c)');

% Plot: Autocorrelation
a_c = xcorr(c_t);
figure(4)
plot(abs(a_c));
title('Autocorrelation of c(nTs)');
legend('Spectrum Method');
xlabel('t');
ylabel('Ac(t)');

% Plot: Power spectral density
% Theoretical - Eqn.5
f = -fd:fd;
psd_th = (1/pi*fd)*(1./((1-(f/fd).^2).^(1/2)));
% Obtained Power Spectral Density
psd_ob = fft(a_c,length(f));
figure(5)
title('Power Spectral Density: Spectrum method');
plot(f,psd_th,'b');
hold on
plot(f,abs(psd_ob),'r');
hold off
legend('Theoretical','Spectrum Method');
xlabel('f');
ylabel('psd(f)');
