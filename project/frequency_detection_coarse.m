function LTF_sig = frequency_detection_coarse(packet_agc, freq_offset_ppm)
% Based on Handout: supl 7 ofdm synch
Fs = 20e6; % Sample Frequency
Fc = 2.4e9; % Carrier Frequency

%%
% Create Freq Offset signal
freq_offset_hz = freq_offset_ppm/1e6*2.4e9;
freq_offset_sig = exp(2j*pi*freq_offset_hz/Fs*(1:length(packet_agc)));

% Add freq offset to packet signal
packet_fdc = packet_agc.*freq_offset_sig;

%% Input Parameters
sig_in = packet_fdc; % Input signal to Coarse Freq Detector

%% Other initializations and pre-set parameters
STF_start_ind = 301;
STF_end_ind = 300+160;
STF_cut = 50;
LTF_len = 64;
% freq_est_hist = zeros(STF_end_ind - 16 ,1);

%% Freq Estimator 

delta_phase = phase(sig_in(STF_start_ind+STF_cut+16:STF_end_ind)*sig_in(STF_start_ind+STF_cut:(STF_end_ind-16))');

figure()
a = sig_in(STF_start_ind+STF_cut+16:STF_end_ind).*sig_in(STF_start_ind+STF_cut:(STF_end_ind-16));
subplot(2,1,1)
plot(abs(a),'.-')
subplot(2,1,2)
plot(mod(phase(a),2*pi),'.-')
title('Debug: FDC');

coarse_freq_offset_est_hz = delta_phase/(2*pi*16/Fs);

coarse_freq_offset_est_ppm = coarse_freq_offset_est_hz/Fc*1e6;
disp(['Freq Offset Actual (ppm): ', num2str(freq_offset_ppm)]);
disp(['Course Freq Offset Est (ppm): ', num2str(coarse_freq_offset_est_ppm)]);

% coarse_freq_offset_sig = exp(-2j*pi*coarse_freq_offset_est_hz/Fs*(1:LTF_len*2));
% LTF_sig = sig_in(STF_end_ind+1:STF_end_ind+LTF_len*2);
% LTF_sig = LTF_sig.*coarse_freq_offset_sig;

coarse_freq_offset_sig = exp(-2j*pi*coarse_freq_offset_est_hz/Fs*(1:numel(sig_in)));
LTF_sig = sig_in.*coarse_freq_offset_sig;
