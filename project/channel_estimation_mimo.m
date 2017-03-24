%% Channel Estimation
% Based on Handout: supl 7 ofdm synch
clear all
close all

d = 100;
SNR_dB = 1000.0;
amplitudes = [0, 0];
delays = [0, 50e-9];

M = 4;

Fs = 20e6; % Sample Frequency
Fc = 2.4e9; % Carrier Frequency
iss = 1;
nss = 1;
n_rx = nss;
tx_packets = create_tx_packet([0],nss,M);

packet11 = create_channel_model(tx_packets{1},amplitudes,delays,d,SNR_dB);

if nss >= 2
packet12 = create_channel_model(tx_packets{2},amplitudes,delays,d,SNR_dB);
packet21 = create_channel_model(tx_packets{1},amplitudes,delays,d,SNR_dB);
packet22 = create_channel_model(tx_packets{2},amplitudes,delays,d,SNR_dB);
end

rx_packets = cell(nss,1);
if nss == 1
rx_packets{1} = packet11;
elseif nss == 2
rx_packets{1} = packet11+packet12;
rx_packets{2} = packet21+packet22;
end

%% Input Parameters

% Input signal is rx_packets

noise_offset = 299;

% LTF1_start_ind = 161; % start of LTF1
LTF1_start_ind = noise_offset+161; % start of LTF1

%% Additional useful variables

GI_len = 16; % Guard Interval length in samples

fft_len = 64;
sym_len = fft_len + GI_len;
NSD_ht = 52;                % NSD for HT
Ntone_htdata = 56;

MIMO_LTF_start_ind = LTF1_start_ind + sym_len*2; % LTF2 starts after LTF1

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

% Order in which data is extracted from matlab fft vector
% [0,1,2,3,...,28,_29,_30,_31,_32,_-31,_30,_-29,-28,-27,...-3,-2,-1]
data_inds = [65-28:65-22,65-20:65-8,65-6:65-1,2:7,9:21,23:29];
data_ks = [-28:-22,-20:-8,-6:-1,1:6,8:20,22:28];

%% Channel Estimator (Only for spatial streams 1-3, spacial stream 4 needs to be negated)

% Source:
% IEEE 802.11n: On Performance of Channel
% Estimation Schemes over OFDM MIMO SpatiallyCorrelated
% Frequency Selective Fading TGn Channels
% Roger Pierre, Fabris Hoefel 

% Initialize Y, H_est, and P
Y = zeros(n_rx,4,NSD_ht);
H_est = zeros(n_rx,nss,NSD_ht);
H_inv_est = zeros(nss,n_rx,NSD_ht);
P = P_htltf(1:nss,:);

% Fill up y matrix
for i_rx = 1:n_rx
    sig = rx_packets{i_rx};
    
    % Extract the two LTF1 symbols
    LTF1_sym1 = sig(LTF1_start_ind + 2*GI_len           : LTF1_start_ind + 2*GI_len +   fft_len - 1);
    LTF1_sym2 = sig(LTF1_start_ind + 2*GI_len + fft_len : LTF1_start_ind + 2*GI_len + 2*fft_len - 1);
    
    LTF2_sym = sig(MIMO_LTF_start_ind +             GI_len : MIMO_LTF_start_ind +             GI_len + fft_len - 1);
    LTF3_sym = sig(MIMO_LTF_start_ind +   sym_len + GI_len : MIMO_LTF_start_ind +   sym_len + GI_len + fft_len - 1);
    LTF4_sym = sig(MIMO_LTF_start_ind + 2*sym_len + GI_len : MIMO_LTF_start_ind + 2*sym_len + GI_len + fft_len - 1);
    
    LTF1_sym1_fft = fft(LTF1_sym1)/fft_len;
    LTF1_sym2_fft = fft(LTF1_sym2)/fft_len;
    
    LTF1_sym_fft = 0.5*(LTF1_sym1_fft + LTF1_sym2_fft);
    LTF2_sym_fft = fft(LTF2_sym)/fft_len;
    LTF3_sym_fft = fft(LTF3_sym)/fft_len;
    LTF4_sym_fft = fft(LTF4_sym)/fft_len;
    
    Y(i_rx,1,:) = LTF1_sym_fft(data_inds);
    Y(i_rx,2,:) = LTF2_sym_fft(data_inds);
    Y(i_rx,3,:) = LTF3_sym_fft(data_inds);
    Y(i_rx,4,:) = LTF4_sym_fft(data_inds);
end

% Fill up H_est matrix and H_inv_est
% Only look at NSD_ht subcarriers
for k = 1:NSD_ht
    H_est(:,:,k) = Y(:,:,k)*P.'/4/expected_LTF(data_inds(k));
    H_inv_est(:,:,k) = inv(H_est(:,:,k));
end

% H_est(H_est==0) = 1;  % handle divide by zeros

% DEPRECATED. Use H_est instead
% Based on Handout: supl 7 ofdm synch, page 34
% Only works with 1 spatial stream
channel = LTF1_sym_fft.*conj(expected_LTF)*56; % magnitude of X_k is 1/56
% channel(channel==0) = 1; % handle divide by zeros


%%
LTF1_sym1_fft_corrected = LTF1_sym1_fft./channel;
LTF1_sym2_fft_corrected = LTF1_sym2_fft./channel;
figure
subplot(2,1,1)
hold all
plot(real(LTF1_sym1_fft),'b.-');
plot(imag(LTF1_sym1_fft),'r.-');
subplot(2,1,2)
hold all
plot(real(LTF1_sym1_fft_corrected),'b.-');
plot(imag(LTF1_sym1_fft_corrected),'r.-');

h11 = H_est(1,1,:);
h11 = h11(:);

figure
% subplot(2,1,1)
hold all
plot(data_ks,real(h11),'.-')
plot(data_ks,imag(h11),'.-')
title('Channel Model')
xlabel('K indices of Freqs')
ylabel('Magnitude');
legend('Real','Imag');
% subplot(2,1,2)
% hold all
% plot(data_ks,real(channel(data_inds)),'.-')
% plot(data_ks,imag(channel(data_inds)),'.-')

% figure
% subplot(2,1,1)
% hold all
% plot(real(LTF1_sym1),'o-')
% plot(real(LTF1_sym2),'o-')
% plot(-real(LTF2_sym),'.-')
% plot(real(LTF3_sym),'.-')
% plot(real(LTF4_sym),'.-')
% subplot(2,1,2)
% hold all
% plot(imag(LTF1_sym1),'o-')
% plot(imag(LTF1_sym2),'o-')
% plot(-imag(LTF2_sym),'.-')
% plot(imag(LTF3_sym),'.-')
% plot(imag(LTF4_sym),'.-')