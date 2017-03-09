close all
clear
% sinc reconstruction

Ts = 0.1;
Fs = 1/Ts;
n = 0:1:40/Ts;
nTs = n*Ts;
y = sin(n*Ts);

% Filter is 201 taps total in this config
numTaps = 100;
sincTaps = -numTaps:numTaps+1;

% tau_t = 0; % Tau_t is in time
% tau_s = tau_t/Ts; % Tau_s is in samples
tau_s = 10;
yprime = sinc_reconstruction(y, tau_s*Ts, Ts);

% tau = 1.5;
% yprime = sinc_interp(n',y,n'+tau);

figure(1);
plot(n,y,'r-o');
hold on
plot(n,yprime,'b-*');
ylim([-11,11])
hold off

% figure(2)
% plot(sinc(tau + sincTaps),'.-')



