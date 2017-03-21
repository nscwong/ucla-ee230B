close all

d = 100;
SNR_dB = 20.0;
amplitudes = [0, 0];
delays = [0, 50e-9];

iss = 1;
nss = 1;
tx_packet = create_tx_packet([],iss,nss);
packet = create_channel_model(tx_packet,amplitudes,delays,d,SNR_dB);

% Create Freq Offset signal
freq_offset_ppm = 50; 
freq_offset_hz = freq_offset_ppm/1e6*2.4e9;
freq_offset_sig = exp(2j*pi*freq_offset_hz/20e6*(1:length(packet)));

% Add freq offset to packet signal
packet = packet.*freq_offset_sig;

%% AGC added below

%% Input Parameters

sig_in = packet; % Input signal to Coarse Freq Detector


%% Other initializations and pre-set parameters
STF_end_ind = 160;

freq_est_hist = zeros(STF_end_ind - 16 ,1);

%% Freq Estimator 

%loop index
ind = 2; % starting at 2 because first sample is weird
% Coarse Freq Detection Loop
while true
   if ind >= STF_end_ind - 16 
       % STF done, STF period length is 16
       break
   end
   
   % Phase difference between 16 sample STF cycles
   delta_phase = phase(conj(sig_in(ind))*sig_in(ind+16));
   % Estimate freq using this
   delta_f_est = delta_phase/(2*pi*0.8e-6);
   
   % store to history
   freq_est_hist(ind) = delta_f_est;  
   
   % Update loop index
   ind = ind + 1;
end
%%
figure
plot(freq_est_hist(2:end)/2.4e9*1e6,'.-')
title('Freq Estimations')
xlabel('Sample')
ylabel('ppm')

coarse_freq_offset_est_ppm = mean(freq_est_hist(2:end))/2.4e9*1e6;
disp(['Course Freq Offset Est (ppm): ', num2str(coarse_freq_offset_est_ppm)];