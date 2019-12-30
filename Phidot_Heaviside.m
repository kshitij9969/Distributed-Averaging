function tdot=Phidot_Heaviside(x)
N=4;
% W=linspace(0.2,0.3,N);
W=0.55;
k=0.27*(ones(N)-eye(N));
d=[0;pi/2;-pi/2;-pi/2];
% d=0;
% d=0;
% m=1;
% for i=1:N-1
%     b=2*pi*i*m/N;
%     d=[d;b];
% end
tdot= W + k*abs(sum(ones(N,1)*heaviside(sin(x))'-heaviside(sin(x+d))*ones(N,1)',2));
end