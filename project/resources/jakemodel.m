% Jake's Model

bandwidth = 20e6;
bandwidth = 1e3;
dt = 1/bandwidth;
t = 0:dt:1; % Time Samples
Fdoppler = 10;

M = 10000; % Number of scatterers
N = 2;   % Number of waveforms
a = unifrnd(0,2*pi,M,N);
b = unifrnd(0,2*pi,M,N);
alpha = unifrnd(0,2*pi,M,N);

R_I = zeros(N, numel(t));
R_Q = zeros(N, numel(t));

for n = 1:N
   R_I(n,:) = sum(cos(2*pi*Fdoppler*cos(alpha(:,n)*t) + a(:,n)*ones(1,numel(t))))/sqrt(M);
   R_Q(n,:) = sum(sin(2*pi*Fdoppler*cos(alpha(:,n)*t) + a(:,n)*ones(1,numel(t))))/sqrt(M);
end

R = R_I + R_Q;

theta = unifrnd(0,2*pi,M,N);


figure;
plot(xcorr(R(1,:)));