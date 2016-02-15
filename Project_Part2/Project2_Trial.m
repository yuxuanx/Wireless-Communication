clc;clear
%% Parameters
N = 32; % number of subcarriers used in OFDM
Nsym = 50; % total number of OFDM symbols
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
s = sqrt(E/2)*(b_buffer(:,1) + b_buffer(:,2)*1j); % bits to symbols
Ncp = 8; % length of cyclic prefix
% vectorize symbols into Nsym blocks, each block contains N symbols
ss = reshape(s,N,Nsym);
z = sqrt(N/Ts)*ifft(ss); % Generate OFDM Seuqence
zCyclic = [z(end-Ncp+1:end,:);z]; % add cyclic prefix
% concatenate OFDM symbols to get transmitted samples
zz = reshape(zCyclic,(N+Ncp)*Nsym,1);
% n = (0:length(zz)-1)';
% ofdmModu = zz.*cos(2*pi*fc/fs*n); % modulation

%% Channel
tau = [0 4]; % delays of taps in samples
fd = v*fc/c; % Doppler frequency
fdTs = fd*Ts; % normalized Doppler frequency
P = [0.5 0.5]; % power delay profile
[r, h] = Fading_Channel(zz, tau, fdTs, P); % generate the fading channel
C = diag(fft(h(1),N));
r = r(1:end-tau(end)); % discard delay samples

%% Receiver
% ofdmDemodu = r.*cos(2*pi*fc/fs*n); % demodulation
y = reshape(r,N+Ncp,Nsym); % serial to parallel
y = y(Ncp+1:end,:); % remove cyclic prefix
rx = sqrt(Ts/N)*fft(y); % FFT 
x = sqrt(N0/2)*(randn(N,Nsym) + 1i*randn(N,Nsym)); % AWGN Channel
bitReceive = zeros(2*M,1); % bits received
sr = C'*rx + x; % symbols received
% decision making, symbols to bits
bitReceive(1:2:end) = reshape(sign(real(sr)),M,1);
bitReceive(2:2:end) = reshape(sign(imag(sr)),M,1);
% calculate the bit error rate
errorRate = length(find((bitReceive - b)~=0))/(2*M);
