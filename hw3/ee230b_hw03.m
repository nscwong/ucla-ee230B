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

N = 1e5; % Number of samples to generate
Beta_stdev = 2; % Standard deviation of log-normal Beta

% Generate log normal vector
Beta_dB = normrnd(0, Beta_stdev, N, 1);
% Convert to linear
Beta = 10.^(Beta_dB/10);

% Plot histograms in linear and log domains
figure(1) 
subplot(2,1,1);
histogram(Beta);
title('Beta')
xlabel('Amplitude')
subplot(2,1,2);
histogram(10*log10(Beta));
title('Beta (dB)')
xlabel('Amplitude (dB)')

% % This is a Rayleigh Distribution (We'll probably need it later)
% C = normrnd(0, 1, T_len, 1) + 1i.*normrnd(0, 1, T_len, 1); % Is this log-normal?
% C_abs = abs(C);
% C_ang = atan2(imag(C), real(C));
% 
% figure; % PDF of amplitude
% histogram(C_abs);
% 
% figure; % PDF of phase
% histogram(C_ang, (-pi):(pi/48):(pi));

%% Part B

% Transmitter uses a square-root-raise-cosine pulse shape
%             generates BPSK modulated signals at rate of 40 Mbps
%             with an average transmit power of 10 dBm
% Center frequency is 2.4 GHz
% Noise PSD is -170 dBm/Hz
% Assume only path loss (no shadowing) between transmitter and receiver
% Receiver uses a sqrt-raised-cosine matched filter to the transmitter

f_c = 2.4 * (10^9); % Center Frequency
alpha = 3; % Path Loss exponent



% Measure and plot the receive SNR after the RX matched filter for TX-RX
% separation distances of 60/80/100/120/140/160 m



% Also plot expected theoretical SNR and explain any differences



% Measure the BER at the same distances as above compare to theory



%% Part C

% Add log normal shadowing to the system of Part B with STDDEV of 4 dB



% Application: requires SNR of 6 dB
% Measure outage probability at a distance of 60/100/140 m
% Compare simulation results with theoretical predictions

