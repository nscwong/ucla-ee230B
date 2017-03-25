% Based on Handout: supl 7 ofdm synch
% Now also includes fine frequency detection.
close all

d = 50;
SNR_dB = 15.0;
amplitudes = [0, 0];
delays = [0, 50e-9];

Fs = 20e6; % Sample Frequency
Fc = 2.4e9; % Carrier Frequency
iss = 1;
nss = 1;
tx_packet = create_tx_packet([0],nss);
packet_agc = create_channel_model(tx_packet{iss},amplitudes,delays,d,SNR_dB);
%%
% Create Freq Offset signal
freq_offset_ppm = 50; 
freq_offset_hz = freq_offset_ppm/1e6*2.4e9;
freq_offset_sig = exp(2j*pi*freq_offset_hz/Fs*(1:length(packet_agc)));

% Add freq offset to packet signal
packet_fdc = packet_agc.*freq_offset_sig;

%% Input Parameters

sig_in = packet_fdc; % Input signal to Coarse Freq Detector


%% Other initializations and pre-set parameters
STF_start_ind = 1;
STF_end_ind = 160;
STF_cut = 20;
LTF_len = 64;
% freq_est_hist = zeros(STF_end_ind - 16 ,1);

%% Freq Estimator 

delta_phase = phase(sig_in(STF_start_ind+STF_cut+16:STF_end_ind)*sig_in(STF_start_ind+STF_cut:(STF_end_ind-16))');
coarse_freq_offset_est_hz = delta_phase/(2*pi*16/Fs);

coarse_freq_offset_est_ppm = coarse_freq_offset_est_hz/Fc*1e6;
disp(['Freq Offset Actual (ppm): ', num2str(freq_offset_ppm)]);
disp(['Course Freq Offset Est (ppm): ', num2str(coarse_freq_offset_est_ppm)]);


coarse_freq_offset_sig = exp(-2j*pi*coarse_freq_offset_est_hz/Fs*(1:LTF_len*2));
LTF_sig = sig_in(STF_end_ind+1:STF_end_ind+LTF_len*2);
LTF_sig = LTF_sig.*coarse_freq_offset_sig;

delta_phase = phase(LTF_sig(LTF_len+1:LTF_len*2)*LTF_sig(1:LTF_len)');
fine_freq_offset_est_hz = delta_phase/(2*pi*LTF_len/Fs);
fine_freq_offset_est_ppm = fine_freq_offset_est_hz/Fc*1e6;

disp(['Fine Freq Offset Est (ppm): ', num2str(fine_freq_offset_est_ppm)]);

total_freq_offset_est_hz = coarse_freq_offset_est_hz+fine_freq_offset_est_hz;
total_freq_offset_est_ppm = total_freq_offset_est_hz/Fc*1e6;

disp(['Total Freq Offset Est (ppm): ', num2str(total_freq_offset_est_ppm)]);

% Total frequency correction
freq_offset_est_sig = exp(-2j*pi*total_freq_offset_est_hz/Fs*(1:length(packet_agc)));
sig_out = sig_in.*freq_offset_est_sig;

% %loop index
% ind = 2; % starting at 2 because first sample is weird
% % Coarse Freq Detection Loop
% while true
%    if ind >= STF_end_ind - 16 
%        % STF done, STF period length is 16
%        break
%    end
%    
%    % Phase difference between 16 sample STF cycles
%    delta_phase = conj(sig_in(ind))*sig_in(ind+16);
%    % TODO Perhaps weight the lower magnitudes less?
%    
%    % Estimate freq using this
%    delta_f_est = delta_phase/(2*pi*0.8e-6);
%    
%    % store to history
%    freq_est_hist(ind) = delta_f_est;  
%    
%    % Update loop index
%    ind = ind + 1;
% end
% %%
% figure
% subplot(2,1,1)
% plot(freq_est_hist(2:end)/2.4e9*1e6,'.-')
% title('Freq Estimations')
% xlabel('Sample')
% ylabel('ppm')
% subplot(2,1,2)
% plot(abs(sig_in(2:length(freq_est_hist)+1)),'.-')% IN PROGRESS

