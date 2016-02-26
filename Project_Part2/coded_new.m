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
EbN0 = 0:1:25; % dB
Ts = 1/fs; % symbol period (minimum according to Sampling Theorem)
%E = Pt*Ts; % energy per transmitted symbol
Eb = N0*10.^(EbN0/10); % energy per bit
E = 2*Eb*pathLoss*N/(N+Ncp);

%Ts = E/Pt;
R = N*m/(N+Ncp)/Ts; % data rate

errorRate = zeros(length(E),1);

iter_num = 1e3; % Monte Carlo
for div = 2:3
    for i = 10:length(E)
        er = zeros(iter_num,1);
        for j = 1:iter_num
            %% Transmitter
            b = randsrc(2*M,1,[-1 1]); % generate data bits to transmit
            % QPSK constellation with Gray-mapping
            b_buffer = buffer(b, m)'; % Group bits into bits per symbol
            s = sqrt(E(i)/2)*(b_buffer(:,1) + b_buffer(:,2)*1j); % bits to symbols
            % vectorize symbols into Nsym blocks, each block contains N symbols
            ss = reshape(s,N,Nsym);
            if div == 2
                ss = [ss ss];
            elseif div == 3
                ss = [ss ss ss];
            end
            z = sqrt(N/Ts)*ifft(ss); % Generate OFDM Seuqence
            zCyclic = [z(end-Ncp+1:end,:);z]; % add cyclic prefix
            % concatenate OFDM symbols to get transmitted samples
            zz = reshape(zCyclic,(N+Ncp)*Nsym*div,1);
            
            sym = zeros(N, Nsym*div);
            for nn = 1:Nsym*div
            %% Channel
            tau = [0 4]; % delays of taps in samples
            fd = v*fc/c; % Doppler frequency
            fdTs = fd*Ts; % normalized Doppler frequency
            % (N+Ncp)*fdTs << 1
            P = [0.5 0.5]; % power delay profile
            [r, h] = Fading_Channel(zCyclic(:, nn), tau, fdTs, P); % generate the fading channel
            hh = [h(1,1) 0 0 0 h(1,2)]; % simulate fading channel gain
            C = diag(fft(hh,N));
            r = r(1:end-tau(end))/sqrt(pathLoss); % discard delay samples
            %% Receiver
%             y = reshape(r,(N+Ncp),Nsym*div); % serial to parallel
            y = r(Ncp+1:end,:); % remove cyclic prefix
            rx = sqrt(Ts/N)*fft(y); % FFT
            x = sqrt(N0/2)*(randn(N,1) + 1i*randn(N,1)); % AWGN Channel
            bitReceive = zeros(2*M,1); % bits received
            sr = C'*(rx + x); % symbols received
            
            sym(:,nn) = sr./(abs(diag(C)).^2);
            sym(:,nn) = sym(:,nn)/sqrt(Eb(i));
%             sym = zeros(N,Nsym*div);
%             for k = 1:Nsym*div
%                 sym(:,k) = sr(:,k)./(abs(diag(C)).^2);
%             end
%             sym = sym/sqrt(Eb(i));
            end
            sr = zeros(N, Nsym);
            if div == 2
                sr = (sym(:,1:Nsym) + sym(:,Nsym+1:Nsym*div))/div;
            elseif div == 3
                sr = (sym(:,1:Nsym) + sym(:,Nsym+1:Nsym*(div-1)) + sym(:,Nsym*(div-1)+1:Nsym*div))/div;
            end
            
            sr = reshape(sr, M, 1);
%             sr = reshape(sym,N*Nsym,div);
%             
%             sr = mean(sr(:,:),2);
            
            Receive = zeros(length(sr),1);
            qpsk = [1+1j,1-1j,-1-1j,-1+1j];
            for kk = 1:length(sr)
                [~,I] = max(real(sr(kk)*conj(qpsk)));
                Receive(kk) = qpsk(I);
            end
            
            bitReceive(1:2:end) = reshape(real(Receive),M,1);
            bitReceive(2:2:end) = reshape(imag(Receive),M,1);
            
            er(j) = length(find((bitReceive - b)~=0));
        end
        errorRate(i,div) = sum(er)/(2*M*iter_num);
    end
end
openfig('codedOFDM.fig');hold on
semilogy(EbN0,errorRate,'.');