%% Parameters
M = 300; % time samples
N = 64; % frequency samples
L = 3; % number of taps

%% 9 combinations for L = 1,2,3 fdts = 0.1,0.01,0.005
channelResp( 0.1, 1, M, N );
channelResp( 0.01, 1, M, N );
channelResp( 0.005, 1, M, N );
channelResp( 0.1, 2, M, N );
channelResp( 0.01, 2, M, N );
channelResp( 0.005, 2, M, N );
channelResp( 0.1, 3, M, N );
channelResp( 0.01, 3, M, N );
channelResp( 0.005, 3, M, N );








