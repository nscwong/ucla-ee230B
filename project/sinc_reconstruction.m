
% sinc reconstruction

Ts = 0.1;
Fs = 1/Ts;
n = 0:1:40/Ts;
nTs = n*Ts;
y = sin(n*Ts);

tau = 1;
yprime = conv(y,sinc(tau + nTs));

figure;
plot(n,y,'-o');
hold on
plot(yprime(1:numel(nTs)),'-*');
hold off