%% Parameters
N = 128; % number of subcarriers used in OFDM
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
b = randsrc(1,2*M,[-1 1]); % generate data bits to transmit
m = 2; % Bits per Symbol in QPSK
% QPSK constellation with Gray-mapping
s_QPSK = sqrt(E/2)*[(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]; 
b_buffer = buffer(b, m)'; % Group bits into bits per symbol
s = b_buffer(:,1) + b_buffer(:,2)*1j;
% initialization
ss = zeros(N,Nsym);
z = zeros(N,Nsym);
Ncp = 10; % length of cyclic prefix
zCyclic = zeros(N+Ncp,Nsym);
for i=1:Nsym
    ss(:,i) = s((i-1)*N+1:i*N);
    z(:,i) = sqrt(N/Ts)*ifft(ss(:,i)); % Generate OFDM Seuqence
    cyclicPrefix = z(end-Ncp+1:end,i); % cyclic prefix
    zCyclic(:,i) = [cyclicPrefix;z(:,i)]; % add cyclic prefix
end
% concatenate OFDM symbols to get transmitted samples
zz = reshape(zCyclic,1,(N+Ncp)*Nsym)';

%% Channel
tau = [0 4]; % delays of taps in samples
fd = v*fc/c; % Doppler frequency
fdTs = fd*Ts; % normalized Doppler frequency
P = [0.5 0.5]; % power delay profile
[r, h] = Fading_Channel(zz, tau, fdTs, P); % generate the fading channel
C = zeros(N,N,Nsym);
for j=1:Nsym
    C(:,:,j) = diag(fft(h((j-1)*N+1,:),N)); % channel response
end

%% Receiver
