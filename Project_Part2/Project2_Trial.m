clc;clear
%% Parameters (Do think twice about the parameters Tb)
N = 64; % number of subcarriers used in OFDM
Nsym = 1; % total number of OFDM symbols
Ncp = 5; % length of cyclic prefix (>=4)
m = 2; % bits per Symbol in QPSK
fc = 2e9; % carrier frequency
fs = 1e6; % sampling frequency
Pt = 0.1; % average transmitted power (W)
pathLoss = 10^(101/10); % linear pass loss
Pr = Pt/pathLoss;
N0 = 2.07e-20*2; % noise spectral density N0/2
v = 15; % speed of receiver
c = 3e8; % speed of light
M = N*Nsym; % number of time samples
EbN0 = 0:1:25; % dB
% Ts = 1/fs; % symbol period (minimum according to Sampling Theorem)
%E = Pt*Ts; % energy per transmitted symbol
Eb = N0*10.^(EbN0/10); % energy per bit
E = 2*Eb*pathLoss*N/(N+Ncp);
ts = E/Pt;
% R = N*m/(N+Ncp)./ts; % data rate

errorRate = zeros(length(E),1);

iter_num = 1e3; % Monte Carlo

parfor i = 1:length(E)
    er = zeros(iter_num,1);
    for j = 1:iter_num
%% Transmitter
b = randsrc(2*M,1,[-1 1]); % generate data bits to transmit
% QPSK constellation with Gray-mapping
b_buffer = buffer(b, m)'; % Group bits into bits per symbol
s = sqrt(E(i)/2)*(b_buffer(:,1) + b_buffer(:,2)*1j); % bits to symbols
% vectorize symbols into Nsym blocks, each block contains N symbols
ss = reshape(s,N,Nsym);
z = sqrt(N/ts(i))*ifft(ss); % Generate OFDM Seuqence
zCyclic = [z(end-Ncp+1:end,:);z]; % add cyclic prefix
% concatenate OFDM symbols to get transmitted samples
zz = reshape(zCyclic,(N+Ncp)*Nsym,1);

%% Channel
Ts = 1/fs;
tau = [0 4]; % delays of taps in samples
fd = v*fc/c; % Doppler frequency
fdTs = fd*ts(i); % normalized Doppler frequency
% (N+Ncp)*fdTs << 1
P = [0.5 0.5]; % power delay profile
[r, h] = Fading_Channel(zz, tau, fdTs, P); % generate the fading channel
hh = [h(1,1) 0 0 0 h(1,2)]; % simulate fading channel gain
C = diag(fft(hh,N));
r = r(1:end-tau(end))/sqrt(pathLoss); % discard delay samples

%% Receiver
y = reshape(r,N+Ncp,Nsym); % serial to parallel
y = y(Ncp+1:end,:); % remove cyclic prefix
rx = sqrt(ts(i)/N)*fft(y); % FFT 
x = sqrt(N0/2)*(randn(N,Nsym) + 1i*randn(N,Nsym)); % AWGN Channel
bitReceive = zeros(2*M,1); % bits received
sr = C'*(rx + x); % symbols received
% decision making, symbols to bits
bitReceive(1:2:end) = reshape(sign(real(sr)),M,1);
bitReceive(2:2:end) = reshape(sign(imag(sr)),M,1);
% calculate the bit error rate
er(j) = length(find((bitReceive - b)~=0));
    end
    errorRate(i) = sum(er)/(2*M*iter_num);
end
figure
semilogy(EbN0,errorRate,'x');hold on
load('ber0.mat');load('ebno0.mat');
semilogy(ebno0,ber0);
xlabel('E_b/N_0 (dB)');ylabel('BER');title('Bit Error Probability')
legend('Empirical','Theory');