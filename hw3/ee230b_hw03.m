% Flat Fading Channel Model
% Narrowband assumes Flat Fading for simplification

% h(t) = Beta(t) * delta(t - tau) where the delay tau is fixed,
% but amplitude Beta(t) is time varying and is often given by a
% colored random process

%% Part A

% Generate Beta(t) assuming white log-normal distribution
% Measure and plot its PDF and compare with theory

t = 1000; % Random value?

C = normrnd(0, 1, t) + 1i.*normrnd(0, 1, t); % Is this log-normal?
C_abs = abs(C);
C_ang = atan2(imag(C), real(C));

figure; % PDF of amplitude
histogram(C_abs, 'Normalization', 'pdf');
title('Flat Fading Channel Model: Amplitude');
xlabel('Amplitude of Complex Gaussian');
ylabel('Probability Density Function');

figure; % PDF of phase
histogram(C_ang, (-pi):(pi/48):(pi), 'Normalization', 'pdf');
title('Flat Fading Channel Model: Phase');
xlabel('Phase of Complex Gaussian (radians)');
ylabel('Probability Density Function');

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
avg_tx_power = 10; % dBm


% Measure and plot the receive SNR after the RX matched filter for TX-RX
% separation distances of 60/80/100/120/140/160 m



% Also plot expected theoretical SNR and explain any differences



% Measure the BER at the same distances as above compare to theory



%% Part C

% Add log normal shadowing to the system of Part B with STDDEV of 4 dB



% Application: requires SNR of 6 dB
% Measure outage probability at a distance of 60/100/140 m
% Compare simulation results with theoretical predictions

