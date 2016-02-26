clc;clear
%% Parameters (Do think twice about the parameters Tb)
N = 64; % number of subcarriers used in OFDM
Nsym = 2; % total number of OFDM symbols
Ncp = 5; % length of cyclic prefix (>=4)
m = 2; % bits per Symbol in QPSK
fc = 2e9; % carrier frequency
fs = 1e6; % sampling frequency
Pt = 0.1; % average transmitted power (W)
pathLoss = 10^(101/10); % linear pass loss
Pr = Pt/pathLoss;
v = 15; % speed of receiver
c = 3e8; % speed of light
M = N*Nsym; % number of time samples
EbN0 = 0:1:20; % dB
Ts = 1/fs; % symbol period (minimum according to Sampling Theorem)

E = Pt*Ts;     %Energy of a symbol
Eb = E/2;      %Energy of a bit

R = N*m/(N+Ncp)/Ts; % data rate
% errorRate = zeros(length(EbN0),1);
errorRate = [];
iter_num = 100; % Monte Carlo

for div = 1:3  %Diversity branches
    for i = 1:length(EbN0)
        er = zeros(iter_num,1);
        N0 = (div*Eb)/10^(EbN0(i)/10); %Calculation of N0
        for j = 1:iter_num
            %% Transmitte
            b = randsrc(2*M,1,[-1 1]); % generate data bits to transmit
            % QPSK constellation with Gray-mapping
            b_buffer = buffer(b, m)'; % Group bits into bits per symbol
            s = sqrt(E/2).*(b_buffer(:,1) + b_buffer(:,2)*1j); % bits to symbols
            % vectorize symbols into Nsym blocks, each block contains N symbols
            ss = reshape(s,N,Nsym);
            z = sqrt(N)*ifft(ss); % Generate OFDM Seuqence
            zCyclic = [z(end-Ncp+1:end,:);z]; % add cyclic prefix
            % concatenate OFDM symbols to get transmitted samples
            zz = reshape(zCyclic,1,(N+Ncp)*Nsym);
            
            %% Channel
            tau = [0 4]; % delays of taps in samples
            fd = v*fc/c; % Doppler frequency
            fdTs = fd*Ts; % normalized Doppler frequency
            % (N+Ncp)*fdTs << 1
            P = [0.5 0.5]; % power delay profile
            
            
            %% Repetition of an OFDM symbol for Diversity
            sr_new = [];
            for l = 1:div
                %Each OFDM symbol passes through the channel
                [r, h] = Fading_Channel(zz, tau, fdTs, P); % generate the fading channel
                hh = [h(1,1) 0 0 0 h(1,2)]; % simulate fading channel gain
                x = sqrt(N0/2).*(randn(length(r),1) + 1i*randn(length(r),1)); % AWGN Channel
                x_n = r + x; % Add AWGN to faded channel
                x_n = x_n(1:end-tau(end)); % discard delay samples
                
                %% Receiver
                x_n = reshape(x_n,N+Ncp,Nsym);
                y = x_n(Ncp+1:end,:); % remove cyclic prefix
                y = reshape(y,N*Nsym,1);
                rx = sqrt(1/N)*fft(y); % FFT
                C = diag(fft(hh,N*Nsym)); % Channel estimation
                sr = C'*rx; % Equalization            
                % Put the OFDM symbols together in a matrix
                sr_new = [sr_new sr];
            end
            
            sr = mean(sr_new(:,:),2);
            for g = 1:length(sr)
                rx1 = sqrt((real(sr(g)-(1+1i))^2) + (imag(sr(g)-(1+1i))^2));
                rx2 = sqrt((real(sr(g)-(1-1i))^2) + (imag(sr(g)-(1-1i))^2));
                rx3 = sqrt((real(sr(g)-(-1+1i))^2) + (imag(sr(g)-(-1+1i))^2));
                rx4 = sqrt((real(sr(g)-(-1-1i))^2) + (imag(sr(g)-(-1-1i))^2));
                [res(g) pos(g)] = min([rx1 rx2 rx3 rx4]);
            end
            
            % decision making, symbols to bits
            array = [1+1i 1-1i -1+1i -1-1i];
            Receive = array(pos);
            
            bitReceive = zeros(2*M,1); % bits received
            % decision making, symbols to bits
            bitReceive(1:2:end) = reshape(real(Receive),M,1);
            bitReceive(2:2:end) = reshape(imag(Receive),M,1);
           
            % calculate the bit error rate
            er(j) = length(find((bitReceive - b)~=0));
        end
        errorRate(i,div) = sum(er)/(2*M*iter_num);
    end
end
openfig('codedOFDM.fig');hold on
semilogy(EbN0,errorRate,'.');