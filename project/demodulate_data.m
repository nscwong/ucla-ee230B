function bits_out = demodulate_data(rx_packets, H_inv_est, LTF_start, n_rx, nss, M, num_bits)

datasym1_start_ind = LTF_start+401-1; % Start after LTF4
GI_len = 16; % Guard Interval length in samples
fft_len = 64; % Number of samples in fft period
sym_len = fft_len+GI_len; % Number of samples in symbol
N_dbps = 52*log2(M); % Number of databits per symbol
NSD_ht = 52;                % NSD for HT
Ntone_htdata = 56;
N_syms = ceil(num_bits/N_dbps); % Number of symbols 

% Order in which data is extracted from matlab fft vector
% [0,1,2,3,...,28,_29,_30,_31,_32,_-31,_30,_-29,-28,-27,...-3,-2,-1]
data_inds = [65-28:65-22,65-20:65-8,65-6:65-1,2:7,9:21,23:29];

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

%     sym_hist = [];

% Data extraction loop
for sym_i = 0:nss:N_syms-1
    stream_Yks = zeros(n_rx,NSD_ht);
    stream_Xks = zeros(nss,NSD_ht);
    % Get symbol from each spatial stream
    for i_rx = 1:n_rx
        % Need to extract all, even if no data

        sig = rx_packets{i_rx}; % select spatial stream signal

        sym_i_adj = floor(sym_i/nss); %symbol in this spatial stream

        % extract time data
        datasym = sig(datasym1_start_ind(i_rx) + sym_i_adj*sym_len + GI_len : datasym1_start_ind(i_rx) + sym_i_adj*sym_len + GI_len +   fft_len - 1);
        % Extract subcarrier components
        datasym_fft = fft(datasym)/fft_len;
        stream_Yks(i_rx,:) = datasym_fft(data_inds); % Extract and store Yks
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

%             sym_hist = horzcat(sym_hist,qam_mod_syms_scaled);

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