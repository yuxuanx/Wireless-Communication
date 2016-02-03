%% Parameters
M = 300; % time samples
N = 64; % frequency samples
L = 3; % number of taps
Ns = 10000;
x = (randn(Ns, 1) + sqrt(-1)*randn(Ns, 1))/sqrt(2); % Gaussian noise
%% 9 combinations for L = 1,2,3 fdts = 0.1,0.01,0.005
channelResp( x, 0.1, 1, M, N );
channelResp( x, 0.01, 1, M, N );
channelResp( x, 0.005, 1, M, N );
channelResp( x, 0.1, 2, M, N );
channelResp( x, 0.01, 2, M, N );
channelResp( x, 0.005, 2, M, N );
channelResp( x, 0.1, 4, M, N );
channelResp( x, 0.01, 7, M, N );
channelResp( x, 0.001 , 3, M, N );








