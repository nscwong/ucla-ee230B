clear;clc;

d = 50;
SNR_dB = 100.0;
% amplitudes = [0];
% delays = [0];
amplitudes = [0, 0];
delays = [0, 50e-9];
%amplitudes = [0, -3, -6, -10];
%delays = [0, 70e-9, 150e-9, 200e-9];
M = 4;

iss = 1;
nss = 1;
tx_packet = create_tx_packet([0],nss, M);
packet = create_channel_model(tx_packet{iss},amplitudes,delays,d,SNR_dB);

%% AGC
% Input Parameters
pow_desired_dB = -6; % Desired power level
P_est_samples = 16; % Number of samples to use to estimate P_avg
AGC_gain_dB = 100; % Initial AGC gain
mu_agc = 0.05; % AGC adaptation parameter
% Other initializations and pre-set parameters
packet_agc = zeros(size(packet)); % Signal output from AGC
Ts = 1/20e6; % Sample time (s)
AGC_off_ind = floor(3e-6/Ts); % Indx at which to turn off AGC, 3 us here
% Arrays for plotting
AGC_gain_history = zeros(size(packet)); % Store AGC gain history
error_history = zeros(size(packet)); % Store error history

%% Packet Detection
threshold = 0.45;
packet_detected = 0;

%% AGC and Packet Detection Loop
time_vec = 1:640;
for t = time_vec
    % -- AGC --
    if 1%t < AGC_off_ind 
        % Apply AGC gain to current sample
        packet_agc(t) = 10^(AGC_gain_dB/20)*packet(t);
        % Estimate power using last P_est_samples 
        pow_est_dB = 10*log10(mean(abs(packet_agc(max(1,t-P_est_samples+1):t)).^2));
        error = pow_est_dB-pow_desired_dB; % error
        % Adjust mu based on how large the error
        if abs(error) < 20
            % We're close, lock in
            mu_agc_adj = mu_agc;
        else
            % We're far, adapt quickly
            mu_agc_adj = 0.2;
        end
        % Update AGC gain using adjusted mu
        AGC_gain_dB = AGC_gain_dB - mu_agc_adj*error;
        % Store history data
        AGC_gain_history(t) = AGC_gain_dB;
        error_history(t) = error;
    else
        % Apply AGC gain to current sample
        packet_agc(t) = 10^(AGC_gain_dB/20)*packet(t);
        AGC_gain_history(t) = AGC_gain_dB;
        error_history(t) = error;
    end
    
    % Saturation of ADC
%     if (real(packet_agc(t)) > 1)
%         packet_agc(t) = 1 + imag(packet_agc(t));
%     elseif real(packet_agc(t)) < -1
%         packet_agc(t) = -1 + imag(packet_agc(t));
%     end
%     if (imag(packet_agc(t)) > 1)
%         packet_agc(t) = 1j + real(packet_agc(t));
%     elseif imag(packet_agc(t)) < -1
%         packet_agc(t) = -1j + real(packet_agc(t));
%     end

    % -- Packet Detection --
    if (t > 80) %&& (packet_detected == 0)
        [packet_detected, ds1, ds2] = packet_detection(packet_agc, iss, nss, t-63, t, threshold);
        if packet_detected && (t <= 300) % Debug mechanism
            fprintf('False Packet Detection at t=%d\n',t);
            packet_detected = 0;
        end
        if packet_detected && (t > 300)
            fprintf('t=%d\n',t); 
            %figure;
            %plot(ds1);
            %title(['t=' num2str(t)]);
        end
    end
end

% -- Frequency Detection Coarse --
freq_offset_ppm = 50; 
LTF_sig = frequency_detection_coarse(packet_agc, freq_offset_ppm);

[test, test2] = htltf_start(packet_agc, iss, nss, min(time_vec), max(time_vec));