function [] = channelResp( x, fdts, L, M, N )
%calculate channel response for discrete-time channel model
C = zeros(N,M);
Cl = zeros(N,M);
for l = 1:L
    tao = (l-1)*2;
    for k = 1:N
        expo = exp(-1j*2*pi*k*tao/N);
        cl = rayleighFading( x, M/2, fdts )';
        Cl(k,:) = cl*expo;
    end
    C = C + Cl;
end
responsePlot( M, N, C );
end

