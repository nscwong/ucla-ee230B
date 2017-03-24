
close all
% 0.25 KByte packet paylaod
bits_per_kB = 1024*8;
num_bits = round(0.25*bits_per_kB);

data_bits = randi([0 1], 1,num_bits);  % 1xn vector of data
M = 16; % M-QAM

nss = 2; % Number of spatial streams
tx_packets = create_tx_packet(data_bits, nss, M);

% sig = tx_packets{1};

datasym1_start_ind = 561; % Start after LTF4
GI_len = 16; % Guard Interval length in samples
fft_len = 64; % Number of samples in fft period
sym_len = fft_len+GI_len; % Number of samples in symbol
N_dbps = 52*log2(M); % Number of databits per symbol

N_syms = ceil(num_bits/N_dbps); % Number of symbols 
%%

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

% Order in which data is extracted from matlab fft vector
% [0,1,2,3,...,28,_29,_30,_31,_32,_-31,_30,_-29,-28,-27,...-3,-2,-1]
inds = [65-28:65-22,65-20:65-8,65-6:65-1,2:7,9:21,23:29];

% Initialize array to store demodulated data
bits_out = zeros(1,N_syms*N_dbps);

% Data extraction loop
for sym_i = 0:nss:N_syms-1
    % Get symbol from each spatial stream
    for iss = 1:nss
        curr_sym = sym_i + (iss -1); % current symbol number
        if curr_sym >= N_syms
            % No more symbols left
            break;
        end
            
        sig = tx_packets{iss}; % select spatial stream signal, currently using transmitted signal
        
        sym_i_adj = floor(sym_i/nss); %symbol in this spatial stream
        
        % extract time data
        datasym = sig(datasym1_start_ind + sym_i_adj*sym_len + GI_len : datasym1_start_ind + sym_i_adj*sym_len + GI_len +   fft_len - 1);
        % Extract subcarrier components
        datasym_fft = fft(datasym)/fft_len;
        % scaled for qam demod
        qam_mod_syms_scaled = datasym_fft(inds)*sqrt(56)/K_mod; 
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

figure
clf
subplot(3,1,1)
hold all
plot(data_bits,'ro-');
plot(bits_out(1:length(data_bits)),'b.-')
ylim([-0.1 1.1]);
subplot(3,1,2)
hold all
plot(data_bits,'ro-');
subplot(3,1,3)
ylim([-0.1 1.1]);
hold all
plot(bits_out(1:length(data_bits)),'b.-')
ylim([-0.1 1.1]);

figure
clf
hold all
plot(real(datasym_fft),'.-');
plot(imag(datasym_fft),'.-');

% plot(real(datasym_fft),'.-')
% plot(imag(datasym_fft),'.-')
%%
figure
clf
hold all
plot(real(sig))
plot(imag(sig))
