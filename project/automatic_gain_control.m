close all
d = 100;
SNR_dB = 20.0;
amplitudes = [0, 0];
delays = [0, 50e-9];

iss = 1;
nss = 1;
tx_packet = create_tx_packet([],iss,nss);
packet = create_channel_model(tx_packet,amplitudes,delays,d,SNR_dB);
t_start = 1;
t_end = 500;

%% AGC added below

%% Input Parameters
pow_desired_dB = -6; % Desired power level

P_est_samples = 16; % Number of samples to use to estimate P_avg
AGC_gain_dB = 100; % Initial AGC gain
mu_agc = 0.05; % AGC adaptation parameter

sig_in = packet; % Input signal to AGC

%% Other initializations and pre-set parameters
sig_out = zeros(size(sig_in)); % Signal output from AGC

Ts = 1/20e6; % Sample time (s)
AGC_off_ind = floor(3e-6/Ts); % Indx at which to turn off AGC, 3 us here

% Arrays for plotting
AGC_gain_history = zeros(size(sig_in)); % Store AGC gain history
error_history = zeros(size(sig_in)); % Store error history

%% AGC
%loop index
ind = 1;
% AGC Loop
while true
   if ind >= AGC_off_ind 
       % Turn off the AGC
       break
   end
   % Apply AGC gain to current sample
   sig_out(ind) = 10^(AGC_gain_dB/20)*sig_in(ind);
   % Estimate power using last P_est_samples 
   pow_est_dB = 10*log10(mean(abs(sig_out(max(1,ind-P_est_samples+1):ind)).^2));
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
   AGC_gain_history(ind) = AGC_gain_dB;
   error_history(ind) = error;
   
   % Update loop index
   ind = ind + 1;
end
% Apply AGC gain to rest of signal
sig_out(ind+1:end) = 10^(AGC_gain_dB/20)*sig_in(ind+1:end);
% Finish filling out history array
AGC_gain_history(ind+1:end) = AGC_gain_dB;
error_history(ind+1:end) = error_history(ind);

%% Plot results
figure
subplot(4,1,1)
plot(AGC_gain_history)
title('AGC Gain History')
subplot(4,1,2)
plot(error_history)
title('AGC Error History')
subplot(4,1,3)
hold all
plot(real(sig_out))
plot(imag(sig_out))
title('Sig_Out')
subplot(4,1,4)
hold all
plot(real(sig_in))
plot(imag(sig_in))
title('Sig_In')