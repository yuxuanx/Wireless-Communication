function [ channelGain ] = rayleighFading( x, N, fdts )
%Rayleigh fading with Clarkes spectrum using filter method
ts = 1e-4; % sample interval(s)
fd = fdts/ts;
%-----impulse response-----%
t = -N*ts:ts:N*ts; % make impulse response causal
Z = 2*pi*fd*abs(t)'; % argument of Bessel function
J = besselj(1/4,Z); % Bessel function of the first kind, order 1/4x
g = J./(abs(t)'.^(1/4));
g(N+1) = (pi*fd)^(1/4)/gamma(5/4);
g_hat = g./sqrt(sum(abs(g).^2)); % normalize impulse response

%-----channel gain-----%
channelGain = conv(x,g_hat); % channel gain
channelGain = channelGain(5001:5300); % just choose 300 samples

end

