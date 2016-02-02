%% Parameters
M = 300; % time samples
N = 64; % frequency samples
L = 3; % number of taps
x = (randn(M, 1) + sqrt(-1)*randn(M, 1))/sqrt(2); % Gaussian noise
%% Tap gains cl(nTs) Rayleigh fading with Clarkes spectrum
len = 10; % the length of window 2*N+1
windowType = 'rectwin';
% window shape/@rectwin @hamming etc.
% channelGain = rayleighFading( x, N, 'hamming', 0.1 );
%% 9 combinations for L = 1,2,3 fdts = 0.1,0.01,0.005
channelResp( x, len, windowType, 0.1, 1, M, N );
channelResp( x, len, windowType, 0.01, 1, M, N );
channelResp( x, len, windowType, 0.005, 1, M, N );
channelResp( x, len, windowType, 0.1, 2, M, N );
channelResp( x, len, windowType, 0.01, 2, M, N );
channelResp( x, len, windowType, 0.005, 2, M, N );
channelResp( x, len, windowType, 0.1, 3, M, N );
channelResp( x, len, windowType, 0.01, 3, M, N );
channelResp( x, len, windowType, 0.005, 3, M, N );








