% Chapter 3 : Example 3.20

%             Reconstruction and aliasing using sinc function

%

% Discrete-time Signal x1(n)

Ts = 0.001; Fs = 1/Ts; n = -10:1:10; nTs = n*Ts;

x = exp(-1000*abs(nTs));
x = sin(nTs*4000);

% Analog Signal reconstruction

Dt = 0.00005;

t = -0.01:Dt:0.01;

xa = x * sinc(Fs*(ones(length(nTs),1)*t-nTs'*ones(1,length(t))));

xb = x * sinc(Fs*(ones(length(nTs),1)*t-nTs'*ones(1,length(t)))) + ...
     x * sinc(Fs*(ones(length(nTs),1)*t-nTs'*ones(1,length(t))));

% check

error = max(abs(xa - exp(-1000*abs(t))));

% Plots

subplot(1,1,1)

subplot(2,1,2);

plot(t*1000,xa);

xlabel('t in msec.'); 
ylabel('xa(t)')

title('Reconstructed Signal from x2(n) using sinc function'); hold on

stem(n*Ts*1000,x); 
plot(t*1000,xb);
hold off