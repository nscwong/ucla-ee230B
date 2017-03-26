% close all
% clear;clc;
%load('rayleigh_channel_ideal.mat');

num_trials = 75;
Fc = 2.4e9; % Carrier Frequency

AWGN_Only_Channel = 0; % Rayleigh if 0
Ideal_AGC = 0;
Ideal_BBD = 0;
Ideal_Channel_Estimation = 1; % No AWGN Noise
Ideal_Frequency_Estimation = 1;

OutputFreqs = 0;
d = 50;
SNR_Values = 0:4:24;
% amplitudes = [0];
% delays = [0];
% amplitudes = [0, 0];
% delays = [0, 50e-9];
amplitudes = [0, -3, -6, -10];
delays = [0, 70e-9, 150e-9, 200e-9];
freq_offset_ppm = 10; %50; 
M = 16;

% Create Filename
if AWGN_Only_Channel
    matname = 'awgn_channel';
else
    matname = 'rayleigh_channel';
end
if Ideal_AGC && Ideal_BBD && Ideal_Channel_Estimation && Ideal_Frequency_Estimation
    matname = [matname '_ideal'];
else
    matname = [matname '_nonideal'];
end
if ~Ideal_AGC
    matname = [matname '_agc'];
end
if ~Ideal_BBD
    matname = [matname '_bbd'];
end
if ~Ideal_Channel_Estimation
    matname = [matname '_chanest'];
end
if ~Ideal_Frequency_Estimation
    matname = [matname '_freqest'];
end
matname = [matname '.mat'];

bits_per_kB = 1024*8;
num_bits = round(0.25*bits_per_kB);

bandwidth = 20e6;
iss = 1;
nss = 1;
n_rx = 2;

if Ideal_Frequency_Estimation
    freq_offset_ppm = 0;
end

% SNR_ErrorRates = zeros(numel(SNR_Values), 1);
% SNR_ErrorCounts = zeros(numel(SNR_Values), num_trials);
% SNR_DroppedPackets = zeros(numel(SNR_Values), 1);
% SNR_MissedHTLTFStarts = zeros(numel(SNR_Values), num_trials);

for snr_idx = 6:6%1:numel(SNR_Values)
SNR_dB = SNR_Values(snr_idx);
total_error_count = 0;
error_count_hist = zeros(num_trials,1);

dropped_packets = 0;
missed_htltf_starts = 0;

for trial = 1:num_trials
disp(['Trial: ', num2str(trial)]);

data_bits = randi([0 1], 1,num_bits);  % 1xn vector of data
tx_packets = create_tx_packet(data_bits,nss,M);

% Create Freq Offset signal
freq_offset_hz = freq_offset_ppm/1e6*2.4e9;
freq_offset_sig = exp(2j*pi*freq_offset_hz/bandwidth*(1:length(tx_packets{1})));

% Apply carrier offset to transmitted signal
for iss = 1:nss
    tx_packets{iss} = tx_packets{iss}.*freq_offset_sig;
end

[packet11, packet11_no_awgn] = create_channel_model(tx_packets{1},amplitudes,delays,d,SNR_dB,AWGN_Only_Channel);
if n_rx == 2
[packet21, packet21_no_awgn] = create_channel_model(tx_packets{1},amplitudes,delays,d,SNR_dB,AWGN_Only_Channel); 
end
if nss >= 2
[packet12, packet12_no_awgn] = create_channel_model(tx_packets{2},amplitudes,delays,d,SNR_dB,AWGN_Only_Channel);
[packet21, packet21_no_awgn] = create_channel_model(tx_packets{1},amplitudes,delays,d,SNR_dB,AWGN_Only_Channel);
[packet22, packet22_no_awgn] = create_channel_model(tx_packets{2},amplitudes,delays,d,SNR_dB,AWGN_Only_Channel);
end

rx_packets = cell(nss,1);
if n_rx == 1 && nss == 1
rx_packets{1} = packet11;
elseif n_rx == 2 && nss == 2
rx_packets{1} = packet11+packet12;
rx_packets{2} = packet21+packet22;
elseif n_rx == 2 && nss == 1
rx_packets{1} = packet11;
rx_packets{2} = packet21;
elseif n_rx == 1 && nss == 2
rx_packets{1} = packet11+packet12;
end

if Ideal_Channel_Estimation
    rx_packets_no_awgn = cell(nss,1);
    if n_rx == 1 && nss == 1
    rx_packets_no_awgn{1} = packet11_no_awgn;
    elseif n_rx == 2 && nss == 2
    rx_packets_no_awgn{1} = packet11_no_awgn+packet12_no_awgn;
    rx_packets_no_awgn{2} = packet21_no_awgn+packet22_no_awgn;
    elseif n_rx == 2 && nss == 1
    rx_packets_no_awgn{1} = packet11_no_awgn;
    rx_packets_no_awgn{2} = packet21_no_awgn;
    elseif n_rx == 1 && nss == 2
    rx_packets_no_awgn{1} = packet11_no_awgn+packet12_no_awgn;
    end
end

start_ind_htltf = zeros(1,n_rx);
corrected_packet = cell(1,n_rx);

for i_rx = 1:n_rx
packet = rx_packets{i_rx};

%% AGC
% Input Parameters
pow_desired_dB = -6; % Desired power level
P_est_samples = 8; % Number of samples to use to estimate P_avg
if Ideal_AGC
    AGC_gain_dB = 180; % Initial AGC gain
else
    AGC_gain_dB = 0; % Initial AGC gain
end
mu_agc = 0.1; % AGC adaptation parameter
% Other initializations and pre-set parameters
packet_agc = zeros(size(packet)); % Signal output from AGC
Ts = 1/20e6; % Sample time (s)
AGC_off_ind = floor(3e-6/Ts); % Indx at which to turn off AGC, 3 us here
% Arrays for plotting
AGC_gain_history = zeros(size(packet)); % Store AGC gain history
error_history = zeros(size(packet)); % Store error history


%% Packet Detection
threshold = 0.45;
if Ideal_BBD
    packet_detected = 301;
else
    packet_detected = 0;
end

%% AGC and Packet Detection Loop
time_vec = 1:numel(rx_packets{1});
if ~Ideal_AGC
for t = time_vec    
    % -- AGC --
    if (~packet_detected) || (packet_detected && (t < AGC_off_ind + packet_detected))
        % Apply AGC gain to current sample
        packet_agc(t) = 10^(AGC_gain_dB/20)*packet(t);
        % Estimate power using last P_est_samples 
        pow_est_dB = 10*log10(mean(abs(packet_agc(max(1,t-P_est_samples+1):t)).^2));
        error = pow_est_dB-pow_desired_dB; % error
        % Adjust mu based on how large the error
        mu_agc_adj = mu_agc;
        
        % Update AGC gain using adjusted mu
        AGC_gain_dB = AGC_gain_dB - mu_agc_adj*error;
        
        if AGC_gain_dB > 250
            AGC_gain_dB = 250;            
        elseif AGC_gain_dB < 0
            AGC_gain_dB = 0;
        end
        
        % Store history data
        AGC_gain_history(t) = AGC_gain_dB;
        error_history(t) = error;
    else
        break;
    end
    
    % Saturation of ADC
    if ~Ideal_AGC
        if (real(packet_agc(t)) > 1)
            packet_agc(t) = 1 + 1j*imag(packet_agc(t));
        elseif real(packet_agc(t)) < -1
            packet_agc(t) = -1 + 1j*imag(packet_agc(t));
        end
        if (imag(packet_agc(t)) > 1)
            packet_agc(t) = 1j + real(packet_agc(t));
        elseif imag(packet_agc(t)) < -1
            packet_agc(t) = -1j + real(packet_agc(t));
        end
    end

    % -- Packet Detection --
    if (t > 80) && (packet_detected == 0)
        [packet_detected, ds1, ds2] = packet_detection(packet_agc, iss, nss, t-63, t, threshold);
        if packet_detected
            packet_detected = t;
        end
    end
end
% Apply AGC-locked gain to rest of packet
packet_agc(t:end) = 10^(AGC_gain_dB/20)*packet(t:end);
AGC_gain_history(t:end) = AGC_gain_dB;
error_history(t:end) = error;
if (~packet_detected) || (packet_detected > 500)
    SNR_DroppedPackets(snr_idx) = SNR_DroppedPackets(snr_idx) + 1;
    total_error_count = total_error_count+num_bits;
    SNR_ErrorCounts(snr_idx, trial) = num_bits;
    SNR_MissedHTLTFStarts(snr_idx, trial) = missed_htltf_starts;
    disp('Packet Error!')
    break;
end
else
    packet_agc = 10^(AGC_gain_dB/20)*packet;
end

if ~Ideal_AGC
% Saturation of ADC
    packet_agc(real(packet_agc) > 1) = 1 + 1j*imag(packet_agc(real(packet_agc) > 1));
    packet_agc(real(packet_agc) < -1) = -1 + 1j*imag(packet_agc(real(packet_agc) < -1));
    packet_agc(imag(packet_agc) > 1) = 1j + real(packet_agc(imag(packet_agc) > 1));
    packet_agc(imag(packet_agc) < -1) = -1j + real(packet_agc(imag(packet_agc) < -1));
end

% -- Frequency Detection Coarse --
if Ideal_Frequency_Estimation
    packet_fc_coarse = packet_agc;
else
    % Add freq offset to packet signal
    if OutputFreqs
        disp(['Freq Offset Actual (ppm): ', num2str(freq_offset_ppm)]);
    end
    coarse_freq_offset_est_hz = frequency_detection_coarse(packet_agc, packet_detected);
    coarse_freq_offset_sig = exp(-2j*pi*coarse_freq_offset_est_hz/bandwidth*(1:numel(packet_agc)));
    if OutputFreqs
        coarse_freq_offset_est_ppm = coarse_freq_offset_est_hz/Fc*1e6;
        disp(['Course Freq Offset Est (ppm): ', num2str(coarse_freq_offset_est_ppm)]);
    end
    packet_fc_coarse = packet_agc.*coarse_freq_offset_sig;
end

% -- HT-LTF Start Index --
if Ideal_BBD
    start_ind_htltf(i_rx) = 460;
else
    [start_ind_htltf(i_rx), ~] = htltf_start(packet_fc_coarse, iss, nss, packet_detected, numel(packet_fc_coarse));
    disp(['HTLTF_start: ', num2str(start_ind_htltf(i_rx))])
    if abs(start_ind_htltf(i_rx) - 460) > 10
        missed_htltf_starts = missed_htltf_starts + 1;
        if start_ind_htltf(i_rx) > 700
           total_error_count = total_error_count+num_bits;
           SNR_ErrorCounts(snr_idx, trial) = num_bits;
           SNR_MissedHTLTFStarts(snr_idx, trial) = missed_htltf_starts;
           disp('HTLTF Error!')
           break;
        end
    end
end

% -- Fine Frequency Detection --
if Ideal_Frequency_Estimation
    packet_fc_fine = packet_fc_coarse;
else
    fine_freq_offset_est_hz = frequency_detection_fine(packet_fc_coarse, start_ind_htltf(i_rx));
    % Total frequency correction
    freq_offset_est_sig = exp(-2j*pi*fine_freq_offset_est_hz/bandwidth*(1:length(packet_agc)));
    if OutputFreqs
        fine_freq_offset_est_ppm = fine_freq_offset_est_hz/Fc*1e6;
        disp(['Fine Freq Offset Est (ppm): ', num2str(fine_freq_offset_est_ppm)]);
    end
    packet_fc_fine = packet_fc_coarse.*freq_offset_est_sig;
    total_freq_offset_est_hz = coarse_freq_offset_est_hz+fine_freq_offset_est_hz;
    if OutputFreqs
        total_freq_offset_est_ppm = total_freq_offset_est_hz/Fc*1e6;
        disp(['Total Freq Offset Est (ppm): ', num2str(total_freq_offset_est_ppm)]);
    end
end

corrected_packet{i_rx} = packet_fc_fine;

end

if (~packet_detected) || (packet_detected > 500) || (start_ind_htltf(i_rx) > 700)
    continue;
end

%%
offset_adj = 8;
if Ideal_Channel_Estimation
    for i = 1:n_rx
        rx_packets_no_awgn{i} = rx_packets_no_awgn{i}*10^(AGC_gain_dB/20);
    end
    [~, H_inv_est] = channel_estimation_final(rx_packets_no_awgn, start_ind_htltf-offset_adj, n_rx, nss);
else
    [~, H_inv_est] = channel_estimation_final(corrected_packet, start_ind_htltf-offset_adj, n_rx, nss);
end

bits_out = demodulate_data(corrected_packet, H_inv_est, start_ind_htltf-offset_adj, n_rx, nss, M, num_bits);

% [~, H_inv_est] = channel_estimation_final(rx_packets, start_ind_htltf, n_rx, nss);
% 
% bits_out = demodulate_data(rx_packets, H_inv_est, start_ind_htltf, n_rx, nss, M, num_bits);

% Output BER stats
[num_errs,ratio] = biterr(data_bits,bits_out(1:length(data_bits)));

% total_error_count = total_error_count+num_errs;
disp(['                                   Error Rate: ', num2str(ratio)])
% error_count_hist(trial) = ratio;

total_error_count = total_error_count+num_errs;
error_count_hist(trial) = ratio;

SNR_ErrorCounts(snr_idx, trial) = num_errs;
SNR_MissedHTLTFStarts(snr_idx, trial) = missed_htltf_starts;

end % trial

total_error_rate = total_error_count/(num_bits*num_trials);
disp(['Total Error Rate: ', num2str(total_error_rate)]);
disp(['Total Error Count: ', num2str(total_error_count)]);
disp(['Dropped Packet Count: ', num2str(dropped_packets)]);
disp(['Missed HTLTF Starts: ', num2str(missed_htltf_starts)]);

%SNR_DroppedPackets(snr_idx) = dropped_packets;
SNR_ErrorRates(snr_idx) = total_error_rate;

end
%%
figure;
clf
semilogy(SNR_Values, SNR_ErrorRates', '.-');
title('BER for AWGN Channel with 1x1 System');
xlabel('SNR Values');
ylabel('Bit Error Rate (BER)');
legend('NonIdeal Simulation: AGC, BBD, CE','Location','best');
%%
save(matname);

%%
% figure
% clf
% subplot(2,1,1)
% hold all
% plot(real(10^(AGC_gain_dB/20)*rx_packets{1}*exp(1j*0)),'ro-')
% plot(real(packet_fc_fine),'b.-')
% % plot(real(corrected_packet{1}),'b.-')
% xlim([900 950])
% subplot(2,1,2)
% hold all
% plot(imag(10^(AGC_gain_dB/20)*rx_packets{1}*exp(1j*0)),'ro-')
% plot(imag(packet_fc_fine),'b.-')
% % plot(imag(corrected_packet{1}),'b.-')
% xlim([900 950])

