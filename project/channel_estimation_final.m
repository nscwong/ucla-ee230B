function [H_est, H_inv_est] = channel_estimation_final(rx_packets, LTF1_start_ind, n_rx, nss)

GI_len = 16; % Guard Interval length in samples

fft_len = 64;
sym_len = fft_len + GI_len;
NSD_ht = 52;                % NSD for HT

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
    LTF1_sym1 = sig(LTF1_start_ind(i_rx) + 2*GI_len           : LTF1_start_ind(i_rx) + 2*GI_len +   fft_len - 1);
    LTF1_sym2 = sig(LTF1_start_ind(i_rx) + 2*GI_len + fft_len : LTF1_start_ind(i_rx) + 2*GI_len + 2*fft_len - 1);
    
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
    H_inv_est(:,:,k) = pinv(H_est(:,:,k));
end