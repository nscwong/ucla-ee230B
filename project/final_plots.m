close all
clear all
d = 100;
SNR_dB = 5000.0;
% amplitudes = [0];
% delays = [0];
amplitudes = [0, 0];
delays = [0, 50e-9];

% 0.25 KByte packet paylaod
bits_per_kB = 1024*8;
num_bits = round(10*0.25*bits_per_kB);
% num_bits = 200;

data_bits = randi([0 1], 1,num_bits);  % 1xn vector of data
M = 64; % M-QAM

nss = 2; % Number of spatial streams
n_rx = nss;
tx_packets = create_tx_packet(data_bits, nss, M);

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

%%
noise_offset = 299;

datasym1_start_ind = noise_offset+561; % Start after LTF4
GI_len = 16; % Guard Interval length in samples
fft_len = 64; % Number of samples in fft period
sym_len = fft_len+GI_len; % Number of samples in symbol
N_dbps = 52*log2(M); % Number of databits per symbol
NSD_ht = 52;                % NSD for HT
Ntone_htdata = 56;
N_syms = ceil(num_bits/N_dbps); % Number of symbols 

%% CHANNEL ESTIMATION

LTF1_start_ind = noise_offset+161; % start of LTF1
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


%% DEMODULATOR



% Scaling factor for data
if M == 4
    K_mod = 1/sqrt(2);
elseif M==16
    K_mod = 1/sqrt(10);
elseif M==64
    K_mod = 1/sqrt(42);
else
    K_mod = 1;
end

%%

% Initialize array to store demodulated data
bits_out = zeros(1,N_syms*N_dbps);

sym_hist = []

% Data extraction loop
for sym_i = 0:nss:N_syms-1
    stream_Yks = zeros(nss,NSD_ht);
    stream_Xks = zeros(nss,NSD_ht);
    % Get symbol from each spatial stream
    for iss = 1:nss
        % Need to extract all, even if no data
        
        sig = rx_packets{iss}; % select spatial stream signal
        
        sym_i_adj = floor(sym_i/nss); %symbol in this spatial stream
        
        % extract time data
        datasym = sig(datasym1_start_ind + sym_i_adj*sym_len + GI_len : datasym1_start_ind + sym_i_adj*sym_len + GI_len +   fft_len - 1);
        % Extract subcarrier components
        datasym_fft = fft(datasym)/fft_len;
        stream_Yks(iss,:) = datasym_fft(data_inds); % Extract and store Yks
    end
    
    % Channel Correction
    for k = 1:NSD_ht
        stream_Xks(:,k) = H_inv_est(:,:,k)*stream_Yks(:,k);
    end
    
    % Bit extraction
    for iss = 1:nss        
        curr_sym = sym_i + (iss -1); % current symbol number
        if curr_sym >= N_syms
            % No more symbols left
            break;
        end
        
        % scaled for qam demod
        qam_mod_syms_scaled = stream_Xks(iss,:)*sqrt(Ntone_htdata)/K_mod; 
                
        sym_hist = horzcat(sym_hist,qam_mod_syms_scaled);
        
        % QAM demodulation
        qam_syms = qamdemod(qam_mod_syms_scaled,M);     
        % Correct bit order
        a = fliplr(de2bi(qam_syms));
        % Reshape into vector
        bits = reshape(a.',1,numel(a));
        % Add to bit list
        bits_out(curr_sym*N_dbps+1:(curr_sym+1)*N_dbps) = bits;
    end
end

% Output BER stats
[number,ratio] = biterr(data_bits,bits_out(1:length(data_bits)))


padamt = N_syms*N_dbps - numel(data_bits);    % Number of zeros to pad
data_bits_padded = [data_bits, zeros(1,padamt)];      % Pad to fit modulation scheme
data_i = reshape(data_bits_padded,log2(M),numel(data_bits_padded)/log2(M));    % Change to int
data_i = (2.^((log2(M)-1):-1:0))*data_i;

data_i_out = reshape(bits_out,log2(M),numel(bits_out)/log2(M));    % Change to int
data_i_out = (2.^((log2(M)-1):-1:0))*data_i_out;

% figure
% clf
% subplot(3,1,1)
% hold all
% plot(data_i,'ro-');
% plot(data_i_out(1:length(data_i)),'b.-')
% % ylim([-0.1 1.1]);
% subplot(3,1,2)
% hold all
% plot(data_i,'ro-');
% subplot(3,1,3)
% % ylim([-0.1 1.1]);
% hold all
% plot(data_i_out(1:length(data_i)),'b.-')
% % ylim([-0.1 1.1]);

figure
clf
hold all
plot(real(sym_hist),imag(sym_hist),'o');
% plot(,'.-');

% plot(real(datasym_fft),'.-')
% plot(imag(datasym_fft),'.-')
%%
% figure
% clf
% hold all
% plot(real(sig))
% plot(imag(sig))
