% Authors: Gourav Khadge, Nathan Wong
% 
% Flat Fading Channel Model
% Narrowband assumes Flat Fading for simplification
% 
% h(t) = Beta(t) * delta(t - tau) where the delay tau is fixed,
% but amplitude Beta(t) is time varying and is often given by a
% colored random process

%% Part A

% Generate Beta(t) assuming white log-normal distribution
% Measure and plot its PDF and compare with theory

N = 1e6; % Number of samples to generate
Beta_stdev = 1; % Standard deviation of log-normal Beta

% Generate log normal vector
Beta_dB = normrnd(0, Beta_stdev, N, 1);
% Convert to linear
Beta = 10.^(Beta_dB/10);

% Plot histograms in linear and log domains
figure(1) 
subplot(2,1,1);
histogram(Beta, 'Normalization', 'pdf');
title('Log Normally Distributed Beta')
xlabel('Amplitude')
ylabel('Probability Density Function');
subplot(2,1,2);
histogram(10*log10(Beta), 'Normalization', 'pdf');
title('Log Normally Distributed Beta (dB)')
xlabel('Amplitude (dB)')
ylabel('Probability Density Function');

% % This is a Rayleigh Distribution (We'll probably need it later)
% 
% t = 1000; % Random value?
% 
% C = normrnd(0, 1, t) + 1i.*normrnd(0, 1, t); % Is this log-normal?
% C_abs = abs(C);
% C_ang = atan2(imag(C), real(C));
% 
% figure; % PDF of amplitude
% histogram(C_abs, 'Normalization', 'pdf');
% title('Flat Fading Channel Model: Amplitude');
% xlabel('Amplitude of Complex Gaussian');
% ylabel('Probability Density Function');
% 
% figure; % PDF of phase
% histogram(C_ang, (-pi):(pi/48):(pi), 'Normalization', 'pdf');
% title('Flat Fading Channel Model: Phase');
% xlabel('Phase of Complex Gaussian (radians)');
% ylabel('Probability Density Function');

%% Reference code for RRC TX and RX code
%https://www.mathworks.com/help/comm/examples/raised-cosine-filtering.html?refresh=true#zmw57dd0e6093
Nsym = 6;           % Filter span in symbol durations
rrc_beta = 0.5;         % Roll-off factor
sampsPerSym = 8;    % Upsampling factor

% Parameters
DataL = 20;             % Data length in symbols
R = 40e6;               % Data rate
Fs = R * sampsPerSym;   % Sampling frequency

% Create a local random stream to be used by random number generators for
% repeatability
hStr = RandStream('mt19937ar', 'Seed', 0);

% Generate random data
x = 2*randi(hStr, [0 1], DataL, 1)-1;
% Time vector sampled at symbol rate in milliseconds
time = 1000 * (0: DataL - 1) / R;

% Design raised cosine filter with given order in symbols
rctFilt3 = comm.RaisedCosineTransmitFilter('Shape', 'Square root','RolloffFactor',rrc_beta,'FilterSpanInSymbols',Nsym, 'OutputSamplesPerSymbol', sampsPerSym);

% Upsample and filter.

yc = rctFilt3([x; zeros(Nsym/2,1)]);

% Filter group delay, since raised cosine filter is linear phase and
% symmetric.
fltDelay = Nsym / (2*R);
% Correct for propagation delay by removing filter transients
yc = yc(fltDelay*Fs+1:end);

% Plot data.
figure(1)
stem(time, x, 'kx'); hold on;
% Plot filtered data.
to = 1000 * (0: DataL*sampsPerSym - 1) / Fs;
plot(to, yc, 'm-'); hold off;
% Set axes and labels.
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Transmitted Data', 'Sqrt. Raised Cosine', 'Location', 'southeast')


% Design and normalize filter.
rcrFilt = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          rrc_beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'InputSamplesPerSymbol',  sampsPerSym, ...
  'DecimationFactor',       1);
% Filter at the receiver.
yr = rcrFilt([yc; zeros(Nsym*sampsPerSym/2, 1)]);
% Correct for propagation delay by removing filter transients
yr = yr(fltDelay*Fs+1:end);
% Plot data.
figure(2)
stem(time, x, 'kx'); hold on;
% Plot filtered data.
plot(to, yr, 'b-'); hold off;
% Set axes and labels.
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Transmitted Data', 'Rcv Filter Output',...
    'Location', 'southeast')

%% Part B

% Transmitter uses a square-root-raised-cosine pulse shape
%             generates BPSK modulated signals at rate of 40 Mbps
%             with an average transmit power of 10 dBm
% Center frequency is 2.4 GHz
% Noise PSD is -170 dBm/Hz
% Assume only path loss (no shadowing) between transmitter and receiver
% Receiver uses a sqrt-raised-cosine matched filter to the transmitter

f_c = 2.4 * (10^9); % Center Frequency
alpha = 3; % Path Loss exponent
noise_psd = -170; % dBm/Hz
avg_tx_power_db = 10; % dBm

avg_tx_power = 10^(avg_tx_power_db/10);

Nsym = 6;           % Filter span in symbol durations
rrc_beta = 0.5;         % Roll-off factor
sampsPerSym = 8;    % Upsampling factor

% Parameters
DataL = 20;             % Data length in symbols
R = 40e6;               % Data rate
Fs = R * sampsPerSym;   % Sampling frequency

% Create a local random stream to be used by random number generators for
% repeatability
hStr = RandStream('mt19937ar', 'Seed', 0);

% Generate random data
bits = 2*randi(hStr, [0 1], DataL, 1)-1;
% Time vector sampled at symbol rate in milliseconds
time = 1000 * (0: DataL - 1) / R;

% Design raised cosine filter with given order in symbols
% Design raised cosine filter with given order in symbols
rctFilt3 = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          rrc_beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'OutputSamplesPerSymbol', sampsPerSym);
% Upsample and filter.

TX = rctFilt3([bits; zeros(Nsym/2,1)]);

% Filter group delay, since raised cosine filter is linear phase and
% symmetric.
fltDelay = Nsym / (2*R);
% Correct for propagation delay by removing filter transients
TX = TX(fltDelay*Fs+1:end);
to = 1000 * (0: DataL*sampsPerSym - 1) / Fs;

TX_Power = sum(TX.^2)*Fs/length(TX);

% Normalize power and set to desired power
TX = TX/sqrt(TX_Power)*sqrt(avg_tx_power);

% Check if new power is desired power
% TX_Power = sum(TX.^2)*Fs/length(TX)

carrier = cos(2*pi*Fs*to)';

passband_TX = TX.*carrier;

figure(1)
plot(passband_TX)

% Measure and plot the receive SNR after the RX matched filter for TX-RX
% separation distances of 60/80/100/120/140/160 m



% Also plot expected theoretical SNR and explain any differences



% Measure the BER at the same distances as above compare to theory



%% Part C

% Add log normal shadowing to the system of Part B with STDDEV of 4 dB



% Application: requires SNR of 6 dB
% Measure outage probability at a distance of 60/100/140 m
% Compare simulation results with theoretical predictions

