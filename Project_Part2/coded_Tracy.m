clc;clear
%-----Parameters-----%
N = 64; % number of subcarriers used in OFDM
Nsym = 10; % total number of OFDM symbols
Ncp = 5; % length of cyclic prefix (>=4)
m = 2; % bits per Symbol in QPSK
fc = 2e9; % carrier frequency
fs = 1e6; % sampling frequency
Pt = 0.1; % average transmitted power (W)
pathLoss = 10^(101/10); % linear pass loss
Pr = Pt/pathLoss; % average received power
N0 = 2.07e-20*2; % noise spectral density N0/2
v = 15; % speed of receiver
c = 3e8; % speed of light
M = N*Nsym; % number of time samples
EbN0 = 0:1:20; % dB
Ts = 1/fs; % symbol period (minimum according to Sampling Theorem)
%E = Pt*Ts; % energy per transmitted symbol
Eb = N0*10.^(EbN0/10); % energy per bit
div = 2; % number of diversity
E = 2*Eb*pathLoss*N/(N+Ncp)/div;
R = N*m/(N+Ncp)/Ts; % data rate
errorRate = zeros(length(E),1);
iter_num = 1e2; % Monte Carlo


for i = 1:length(E)
    er = zeros(iter_num,1);
    for j = 1:iter_num
        sr = zeros(N*Nsym,div);
        
        %-----Transmitter-----%
        b = randsrc(2*M,1,[-1 1]); % generate data bits to transmit
        % QPSK constellation with Gray-mapping
        b_buffer = buffer(b, m)'; % Group bits into bits per symbol
        s = sqrt(E(i)/2)*(b_buffer(:,1) + b_buffer(:,2)*1j); % bits to symbols
        % vectorize symbols into Nsym blocks, each block contains N symbols
        ss = reshape(s,N,Nsym);
        
        z = sqrt(N/Ts)*ifft(ss); % Generate OFDM Seuqence
        zCyclic = [z(end-Ncp+1:end,:);z]; % add cyclic prefix
        % concatenate OFDM symbols to get transmitted samples
        zz = reshape(zCyclic,(N+Ncp)*Nsym,1);
        
        %----- Channel-----%
        tau = [0 4]; % delays of taps in samples
        fd = v*fc/c; % Doppler frequency
        fdTs = fd*Ts; % normalized Doppler frequency
        % (N+Ncp)*fdTs << 1
        P = [0.5 0.5]; % power delay profile
        
        % using for loop to achieve time diversity
        for dd = 1:div
            [r, h] = Fading_Channel(zz, tau, fdTs, P); % generate the fading channel
            hh = [h(1,1) 0 0 0 h(1,2)]; % simulate fading channel gain
            C = diag(fft(hh,N));
            r = r(1:end-tau(end))/sqrt(pathLoss); % discard delay samples
            
            %-----Receiver-----%
            y = reshape(r,(N+Ncp),Nsym); % serial to parallel
            y = y(Ncp+1:end,:); % remove cyclic prefix
            rx = sqrt(Ts/N)*fft(y); % FFT
            x = sqrt(N0/2)*(randn(N,Nsym) + 1i*randn(N,Nsym)); % AWGN Channel
            bitReceive = zeros(2*M,1); % bits received
            sr(:,dd) = reshape(C'*(rx+x),M,1); % symbols received
        end
        srr = mean(sr,2)/sqrt(Eb(i));
        Receive = zeros(length(srr),1);
        % QPSK soft demodulator
        qpsk = [1+1j,1-1j,-1-1j,-1+1j];
        for kk = 1:length(srr)
            [~,I] = max(real(srr(kk)*conj(qpsk)));
            Receive(kk) = qpsk(I);
        end
        
        % decision making, symbols to bits
        bitReceive(1:2:end) = real(Receive);
        bitReceive(2:2:end) = imag(Receive);
        er(j) = length(find((bitReceive - b)~=0))/(2*M);
        errorRate(i) = mean(er);
    end
end
openfig('codedOFDM.fig');hold on
semilogy(EbN0,errorRate,'.');