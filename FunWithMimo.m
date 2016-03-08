clear all;
close all;

% note every time you run, you will get a new channel and new figures.
% Sometime the channel will be good, sometime bad

% example of equalization
Mt=2;                       % plotting will only work for Mt = 2
Mr=2;                       % this you can set how you want
rho=50;                     % SNR
H=randn(Mr,Mt);             % the Mr x Mt MIMO channel 

QZF=inv(H'*H)*H';           % ZF equalizer 
QMMSE=inv(H'*H+eye(2)/rho)*H';  % MMSE equalizer 
[U,S,V]=svd(H);             % SVD 

sigma=diag(S);
singularValues=sigma(1:min(Mt,Mr))         % singular values of usable channels


N=1000;                     % number of Monte Carlo runs

for k=1:N
    x=2*round(rand(Mt,1))-1;    % BPSK data
    n=randn(Mr,1);              % AWGN
    
    % for CSIR        
    y(:,k)=sqrt(rho)*H*x+n;
    zZF(:,k)=QZF*y(:,k);
    zMMSE(:,k)=QMMSE*y(:,k);
    
    % for CSIT   
    x_precoded=sqrt(rho)*V*x;
    y_ch=H*x_precoded + n;
    y_shaped=U'*y_ch;
    ySVD(:,k)=y_shaped(1:min(Mt,Mr));   % since only the first Mt streams can be used   
end

% scale the outputs to make them fit in one plot
y=y./norm(y)/N;
zZF=zZF./norm(zZF)/N;
zMMSE=zMMSE./norm(zMMSE)/N;

if (Mt<=2)
    plot(zZF(1,:),zZF(2,:),'g*',zMMSE(1,:),zMMSE(2,:),'b.',y(1,:),y(2,:),'r.')
    axis equal
    title('ZF: green   MMSE: blue   orginal observation: red')
    figure
    plot(ySVD(1,:),ySVD(2,:),'.')
    grid
    axis equal
    title('2 equivalent SVD channels (overlap means that one singular value is very small)')    
end