% Authors: Gourav Khadge, Nathan Wong
% 
% Flat Fading Channel Model
% Narrowband assumes Flat Fading for simplification
% 
% h(t) = Beta(t) * delta(t - tau) where the delay tau is fixed,
% but amplitude Beta(t) is time varying and is often given by a
% colored random process

nbits = 20;	% Number of bits to send
d = 60;        % Distance between TX and RX (m)

%% Part A
% Generate Beta(t) assuming white log-normal distribution
% Measure and plot its PDF and compare with theory

Beta_stdev = 1; % Standard deviation of log-normal Beta

% Generate log normal vector
Beta_dB = normrnd(0, Beta_stdev, nbits, 1);
% Convert to linear
Beta = 10.^(Beta_dB/10);

% Plot histograms in linear and log domains
% figure(1) 
% subplot(2,1,1);
% histogram(Beta, 'Normalization', 'pdf');
% title('Log Normally Distributed Beta')
% xlabel('Amplitude')
% ylabel('Probability Density Function');
% subplot(2,1,2);
% histogram(10*log10(Beta), 'Normalization', 'pdf');
% title('Log Normally Distributed Beta (dB)')
% xlabel('Amplitude (dB)')
% ylabel('Probability Density Function');

%% Part B
% Transmitter uses a square-root-raised-cosine pulse shape
%             generates BPSK modulated signals at rate of 40 Mbps
%             with an average transmit power of 10 dBm
% Center frequency is 2.4 GHz
% Noise PSD is -170 dBm/Hz
% Assume only path loss (no shadowing) between transmitter and receiver
% Receiver uses a sqrt-raised-cosine matched filter to the transmitter

% Measure and plot the receive SNR after the RX matched filter for TX-RX
% separation distances of 60/80/100/120/140/160 m
% Also plot expected theoretical SNR and explain any differences
% Measure the BER at the same distances as above compare to theory

% PARAMETERS
data_rate = 40e6;       % Signal data rate
f_c = 2.4e9;            % Center Carrier Frequency
TX_power_avg_dBm = 10;  % dBm
noise_psd_dBm = -170;   % Noise Power Spectral Density

% GENERATE BITS AND RESPECTIVE SYMBOLS
TX_bits = randi([0,1], nbits, 1);
TX_syms = 2*TX_bits - 1;    % BPSK

% GENERATE SQRT-RAISED-COSINE FILTERED DATA
rrc_nsym = 8;           % Filter span in symbol durations
rrc_beta = 0;           % Roll-off factor
rrc_sampspersym = 4;    % Upsampling factor
Fs = data_rate * rrc_sampspersym;   % Sampling frequency
rrc_bw = Fs*(1+rrc_beta);    % Bandwidth (Fs or data rate?)

filter_group_delay = rrc_nsym/(2*data_rate);
sample_start_point = filter_group_delay*Fs+1;

% Generate filter by MATLAB package
rctFilt = comm.RaisedCosineTransmitFilter( ...
  'Shape',                  'Square root', ...
  'RolloffFactor',          rrc_beta, ...
  'FilterSpanInSymbols',    rrc_nsym, ...
  'OutputSamplesPerSymbol', rrc_sampspersym);
%fvtool(rctFilt, 'Analysis', 'impulse')

% Create samples for TX
TX = step(rctFilt, [TX_syms; zeros(rrc_nsym/2,1)]);
TX_power = sum(TX.^2)/length(TX);
TX_power_avg = 10^(TX_power_avg_dBm/10);
TX = TX/sqrt(TX_power)*sqrt(TX_power_avg);
TX_power_avg_new = sum(TX.^2)/length(TX);

TX_samples = TX(sample_start_point:end);
TX_samples = TX_samples(1:rrc_sampspersym:end);

% figure;
% stem(1:numel(TX), TX);
% hold on
% stem(sample_start_point:rrc_sampspersym:numel(TX), TX_samples);
% hold off

% GENERATE CHANNEL

% Part B: Path Loss only model
alpha = 3;          % Path loss exponent
lambda = 3e8/f_c;   % wavelength, c = 3e8 m/s speed of light
path_gain = (lambda/(4*pi))^2 * (1/d)^alpha;
RX_B_nonoise = TX*path_gain;

% Generate noise vector
noise_psd_dB = noise_psd_dBm - 30;  % dB/Hz
noise_psd = 10^(noise_psd_dB/10);   % Watts/Hz
noise_power = rrc_bw*noise_psd;
noise_power_mW = noise_power*10^3;
noise_power_dB = 10*log10(noise_power);
noise_power_dBm = noise_power_dB + 30;
Noise_1 = randn(size(TX))*sqrt(noise_power_mW);

noise_snr_dB = TX_power_avg_dBm - noise_power_dBm;
noise_snr = 10^(noise_snr_dB/10);
N0 = TX_power_avg_new/noise_snr;
Noise_2 = randn(size(TX))*sqrt(N0);

RX_B = RX_B_nonoise + Noise_1;

figure;
stem(1:numel(RX_B_nonoise), RX_B_nonoise);
hold on
stem(1:numel(RX_B), RX_B);
hold off