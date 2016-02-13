clc;clear
%% Parameters
N = 32; % number of subcarriers used in OFDM
Nsym = 100; % total number of OFDM symbols
fc = 2e9; % carrier frequency
fs = 1e6; % sampling frequency
Pt = 0.1; % average transmitted power
N0 = 2.07e-20*2; % noise spectral density N0/2
Ts = 1/fs; % sample period
E = Pt*Ts; % energy per transmitted symbol
Tb = Ts/2; % time per bit using QPSK
Rb = 1/Tb; % data rate
pathLoss = 10^(101/10); % linear pass loss
Pr = Pt/pathLoss; % average received power
Eb = Pr*Tb; % average received energy per bit
v = 15; % speed of receiver
c = 3e8; % speed of light
M = N*Nsym; % number of time samples

%% Transmitter
b = randsrc(2*M,1,[-1 1]); % generate data bits to transmit
m = 2; % Bits per Symbol in QPSK
% QPSK constellation with Gray-mapping
b_buffer = buffer(b, m)'; % Group bits into bits per symbol
s = sqrt(E/2)*(b_buffer(:,1) + b_buffer(:,2)*1j);
% initialization
ss = zeros(N,Nsym);
z = zeros(N,Nsym);
Ncp = 8; % length of cyclic prefix
zCyclic = zeros(N+Ncp,Nsym);
% vectorize symbols into Nsym blocks, each block contains N symbols
for i=1:Nsym
    ss(:,i) = s((i-1)*N+1:i*N);
    z(:,i) = sqrt(N/Ts)*ifft(ss(:,i)); % Generate OFDM Seuqence
    cyclicPrefix = z(end-Ncp+1:end,i); % cyclic prefix
    zCyclic(:,i) = [cyclicPrefix;z(:,i)]; % add cyclic prefix
end
% concatenate OFDM symbols to get transmitted samples
zz = reshape(zCyclic,1,(N+Ncp)*Nsym)';
% n = (0:length(zz)-1)';
% ofdmModu = real(zz.*exp(1i*2*pi*fc/fs*n)); % modulation

%% Channel
tau = [0 4]; % delays of taps in samples
fd = v*fc/c; % Doppler frequency
fdTs = fd*Ts; % normalized Doppler frequency
P = [0.5 0.5]; % power delay profile
[r, h] = Fading_Channel(zz, tau, fdTs, P); % generate the fading channel
C = diag(fft(h(1)',N));
x = sqrt(N0/2)*(randn(length(r),1) + 1i*randn(length(r),1)); % AWGN Channel
r = r + x;
r = r(1:end-tau(end)); % discard delay samples

%% Receiver
% ofdmDemodu = r.*exp(-1i*2*pi*fc/fs*n); % demodulation
% % anti-aliasing filter
% antiFilter = firpm(32,2*[0 .2 .3 .5],[1 1 0 0]);
% y = conv(ofdmDemodu,antiFilter); % y is now the interpolated signal
% y = y(17:end-16);

rr = reshape(r,N+Ncp,Nsym); % serial to parallel
rr = rr(Ncp+1:end,:); % remove cyclic prefix
rx = sqrt(Ts/N)*fft(conj(rr)); % FFT (note: take conjugate value here)
b_hat = zeros(2*N,Nsym);
% decision making, symbols to bits
for i = 1:Nsym
b_hat(1:2:end,i) = sign(real(C'*rx(:,i)));
b_hat(2:2:end,i) = sign(imag(C'*rx(:,i)));
end
bitReceive = reshape(b_hat,2*N*Nsym,1);
% calculate the bit error rate
errorRate = length(find((bitReceive - b)~=0))/(2*N*Nsym);