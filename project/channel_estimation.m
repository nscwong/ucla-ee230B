%% Channel Estimation
% Based on Handout: supl 7 ofdm synch

close all

d = 100;
SNR_dB = 100.0;
amplitudes = [0, 0];
delays = [0, 50e-9];

Fs = 20e6; % Sample Frequency
Fc = 2.4e9; % Carrier Frequency
iss = 1;
nss = 1;
tx_packet = create_tx_packet([],iss,nss);
packet = create_channel_model(tx_packet,amplitudes,delays,d,SNR_dB);


%% Input Parameters

sig_in = packet; % Input signal to Coarse Freq Detector

LTF_start_ind = 161+32; % start after guard interval
LTF_len = 64;

%% Channel Estimator (Only for spatial streams 1-3, spacial stream 4 needs to be negated)

% Extract the two LTF1 symbols
LTF1_sym1 = sig_in(LTF_start_ind:LTF_start_ind+LTF_len-1);
LTF1_sym2 = sig_in(LTF_start_ind+LTF_len:LTF_start_ind+LTF_len*2-1);

LTF1_sym1_fft = fft(LTF1_sym1);
LTF1_sym2_fft = fft(LTF1_sym2);
% a = fftshift(fft(LTF_sym1));
% b = fftshift(fft(LTF_sym2));

% manually shifted so k =
% [0,1,2,3,...,28,_29,_30,_31,_32,_-31,_30,_-29,-28,-27,...-3,-2,-1]
% to line up with output of fft
expected_LTF1 = [0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1,0,0,0,0,0,0,0,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1];

% Based on Handout: supl 7 ofdm synch, page 34
channel = 0.5*(LTF1_sym1_fft+LTF1_sym2_fft).*conj(expected_LTF1);
channel(channel==0) = 1; % handle divide by zeros


%%
LTF1_sym1_fft_corrected = LTF1_sym1_fft./channel;
LTF1_sym2_fft_corrected = LTF1_sym2_fft./channel;
figure
subplot(2,1,1)
hold all
plot(real(LTF1_sym1_fft),'r.-');
plot(imag(LTF1_sym1_fft),'b.-');
subplot(2,1,2)
hold all
plot(real(LTF1_sym1_fft_corrected),'r.-');
plot(imag(LTF1_sym1_fft_corrected),'b.-');
