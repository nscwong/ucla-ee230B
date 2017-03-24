%% Channel Estimation
% Based on Handout: supl 7 ofdm synch
clear all
close all

d = 100;
SNR_dB = 100.0;
amplitudes = [0, 0];
delays = [0, 50e-9];

Fs = 20e6; % Sample Frequency
Fc = 2.4e9; % Carrier Frequency
iss = 1;
nss = 1;
n_rx = nss;
tx_packets = create_tx_packet([],nss);

packet11 = create_channel_model(tx_packets{1},amplitudes,delays,d,SNR_dB);
% packet12 = create_channel_model(tx_packets{2},amplitudes,delays,d,SNR_dB);
% packet21 = create_channel_model(tx_packets{1},amplitudes,delays,d,SNR_dB);
% packet22 = create_channel_model(tx_packets{2},amplitudes,delays,d,SNR_dB);

rx_packets = cell(nss,1);
% rx_packets{1} = packet11+packet12;
% rx_packets{2} = packet21+packet22;

rx_packets{1} = packet11;
% rx_packets{1} = tx_packets{1};

%% Input Parameters

% Input signal is rx_packets

LTF1_start_ind = 300+161; % start of LTF1

%% Additional useful variables

GI_samples = 16; % Guard Interval length in samples

sym_len = 64;
LTF_len = sym_len + GI_samples;

MIMO_LTF_start_ind = LTF1_start_ind + LTF_len*2; % LTF2 starts after LTF1

% Polarity matrix -- Page 27 and 28 of spec
P_htltf = [ ...     % 
     1 -1  1  1;    % P_htlf(i,n) represents polarity of the
     1  1 -1  1;    % ith spatial stream in the
     1  1  1 -1;    % nth HT training symbol
    -1  1  1  1;    %
];

% manually shifted so k =
% [0,1,2,3,...,28,_29,_30,_31,_32,_-31,_30,_-29,-28,-27,...-3,-2,-1]
% to line up with output of fft
expected_LTF = 1/sqrt(56)*[0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1,0,0,0,0,0,0,0,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1];


%% Channel Estimator (Only for spatial streams 1-3, spacial stream 4 needs to be negated)

% Source:
% IEEE 802.11n: On Performance of Channel
% Estimation Schemes over OFDM MIMO SpatiallyCorrelated
% Frequency Selective Fading TGn Channels
% Roger Pierre, Fabris Hoefel 

% Initialize Y, H_est, and P
Y = zeros(n_rx,4,64);
H_est = zeros(n_rx,nss,64);
P = P_htltf(1:nss,:);

% Fill up y matrix
for i_rx = 1:n_rx
    sig = rx_packets{i_rx};
    
    % Extract the two LTF1 symbols
    LTF1_sym1 = sig(LTF1_start_ind + 2*GI_samples           : LTF1_start_ind + 2*GI_samples +   sym_len - 1);
    LTF1_sym2 = sig(LTF1_start_ind + 2*GI_samples + sym_len : LTF1_start_ind + 2*GI_samples + 2*sym_len - 1);
    
    LTF2_sym = sig(MIMO_LTF_start_ind +             GI_samples : MIMO_LTF_start_ind +             GI_samples + sym_len - 1);
    LTF3_sym = sig(MIMO_LTF_start_ind +   LTF_len + GI_samples : MIMO_LTF_start_ind +   LTF_len + GI_samples + sym_len - 1);
    LTF4_sym = sig(MIMO_LTF_start_ind + 2*LTF_len + GI_samples : MIMO_LTF_start_ind + 2*LTF_len + GI_samples + sym_len - 1);
    
    LTF1_sym1_fft = fft(LTF1_sym1);
    LTF1_sym2_fft = fft(LTF1_sym2);
    
    LTF1_sym_fft = 0.5*(LTF1_sym1_fft + LTF1_sym2_fft);
    LTF2_sym_fft = fft(LTF2_sym);
    LTF3_sym_fft = fft(LTF3_sym);
    LTF4_sym_fft = fft(LTF4_sym);
    
    Y(i_rx,1,:) = LTF1_sym_fft;
    Y(i_rx,2,:) = LTF2_sym_fft;
    Y(i_rx,3,:) = LTF3_sym_fft;
    Y(i_rx,4,:) = LTF4_sym_fft;
end

% Fill up H_est matrix
for k = 1:64
    H_est(:,:,k) = Y(:,:,k)*P.'/4/expected_LTF(k);
end

% H_est(H_est==0) = 1;  % handle divide by zeros

% Based on Handout: supl 7 ofdm synch, page 34
% Only works with 1 spatial stream
channel = LTF1_sym_fft.*conj(expected_LTF);
% channel(channel==0) = 1; % handle divide by zeros


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

h11 = H_est(1,1,:);
h11 = h11(:);

figure
subplot(2,1,1)
hold all
plot(real(h11))
plot(imag(h11))
subplot(2,1,2)
hold all
plot(real(channel))
plot(imag(channel))