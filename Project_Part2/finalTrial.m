clc;clear
%% Parameters (Do think twice about the parameters Tb)
N = 256; % number of subcarriers used in OFDM
Nsym = 100; % total number of OFDM symbols
Ncp = 4; % length of cyclic prefix (>=4)
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
EbN0 = 0:1:20; % dB
Ts = 1/fs; % symbol period (minimum according to Sampling Theorem)
%E = Pt*Ts; % energy per transmitted symbol
Eb = N0*10.^(EbN0/10); % energy per bit
E = 2*Eb*pathLoss*N/(N+Ncp);
% R = N*m/(N+Ncp)./ts; % data rate

%-----Repetition and Interleaving-----%
rp = 1/2; % code rate used in repetition code
index_matrix = reshape(1:2*M/rp,2*M,1/rp);
index = reshape(index_matrix',2*M/rp,1); % index used for interleaving


errorRate = zeros(length(E),1);
iter_num = 1e3; % Monte Carlo

for i = 1:length(E)
    er = zeros(iter_num,1);
    for j = 1:iter_num
%% Transmitter
b = randsrc(2*M,1,[-1 1]); % generate data bits to transmit

bitsRep = reshape(repmat(b,1,1/rp),2*M/rp,1); % repetition encoder
bitsRep(index) = bitsRep; % interleaving

% QPSK constellation with Gray-mapping
b_buffer = buffer(bitsRep, m)'; % Group bits into bits per symbol
s = sqrt(E(i)/2)*(b_buffer(:,1) + b_buffer(:,2)*1j); % bits to symbols
% vectorize symbols into Nsym blocks, each block contains N symbols
ss = reshape(s,N,Nsym/rp);
z = sqrt(N/Ts)*ifft(ss); % Generate OFDM Seuqence
zCyclic = [z(end-Ncp+1:end,:);z]; % add cyclic prefix
% concatenate OFDM symbols to get transmitted samples
zz = reshape(zCyclic,(N+Ncp)*Nsym/rp,1);

%% Channel
% Ts = 1/fs;
tau = [0 4]; % delays of taps in samples
fd = v*fc/c; % Doppler frequency
fdTs = fd*Ts; % normalized Doppler frequency
% (N+Ncp)*fdTs << 1
P = [0.5 0.5]; % power delay profile
[r, h] = Fading_Channel(zz, tau, fdTs, P); % generate the fading channel
hh = [h(1,1) 0 0 0 h(1,2)]; % simulate fading channel gain
C = diag(fft(hh,N));
r = r(1:end-tau(end))/sqrt(pathLoss); % discard delay samples

%% Receiver
y = reshape(r,N+Ncp,Nsym/rp); % serial to parallel
y = y(Ncp+1:end,:); % remove cyclic prefix
rx = sqrt(Ts/N)*fft(y); % FFT 
x = sqrt(N0/2)*(randn(N,Nsym/rp) + 1i*randn(N,Nsym/rp)); % AWGN Channel
sr = C'*(rx+x); % symbols received

%-----soft decision making-----%
% sym = zeros(N,Nsym/rp);
% for k = 1:Nsym/rp
% sym(:,k) = sr(:,k)*sqrt(pathLoss)./(abs(diag(C)).^2);
% end
sym = reshape(sr,M/rp,1)./repmat(abs(diag(C)).^2,Nsym/rp,1);
symNorm = sym/sqrt(Eb(i));

symbols = zeros(2*length(symNorm),1);
symbols(1:2:end) = real(symNorm);
symbols(2:2:end) = imag(symNorm);
symbolsDe = symbols(reshape(index_matrix',1,2*M/rp));
symbolsDe = reshape(symbolsDe, length(b), 1/rp);
symbolsDe = mean(symbolsDe, 2);
symbolsDe = reshape(symbolsDe, 2, length(b)/2);
symbolsDe = symbolsDe(1,:) + 1j*symbolsDe(2,:);

qpsk = [1+1j,1-1j,-1-1j,-1+1j];
symML = zeros(M,1);
for kk = 1:M
    [~,I] = max(real(symbolsDe(:,kk)*conj(qpsk)));
    symML(kk) = qpsk(I);
end
bitReceive = zeros(2*M,1); % bits received
% demap symbols to bits
bitReceive(1:2:end) = sign(real(symML));
bitReceive(2:2:end) = sign(imag(symML));
% deinterleaving
% bitReceive = bitReceive(reshape(index_matrix',1,2*M/rp));
% % maximum ratio combining
% bitDecoded = sum(reshape(bitReceive,2*M,1/rp),2);
% % decode repetition code
% bitDecoded(bitDecoded<0) = -1;
% bitDecoded(bitDecoded>=0) = 1;

er(j) = length(find((bitReceive - b)~=0));
    end
    errorRate(i) = sum(er)/(2*M*iter_num);
end
%% Plotting
openfig('codedOFDM.fig');hold on
semilogy(EbN0,errorRate,'.');