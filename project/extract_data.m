
close all
% 0.25 KByte packet paylaod
bits_per_kB = 1024*8;
num_bits = round(0.25*bits_per_kB);

data_bits = randi([0 1], 1,num_bits);  % 1xn vector of data
M = 64; %QPSK

nss = 1;
tx_packets = create_tx_packet(data_bits, nss, M);

sig = tx_packets{1};

datasym1_start_ind = 561; % Start after LTF4
GI_len = 16; % Guard Interval length in samples
fft_len = 64;
sym_len = fft_len+GI_len;
N_dbps = 52*log2(M);

N_syms = ceil(num_bits/N_dbps);

if M == 4
    K_mod = 1/sqrt(2);
elseif M==16
    K_mod = 1/sqrt(10);
elseif M==64
    K_mod = 1/sqrt(42);
else
    K_mod = 1;
end

% manually shifted so k =
% [0,1,2,3,...,28,_29,_30,_31,_32,_-31,_30,_-29,-28,-27,...-3,-2,-1]
inds = [65-28:65-22,65-20:65-8,65-6:65-1,2:7,9:21,23:29];
% k = [0:32,-31:-1]

bits_out = zeros(1,N_syms*N_dbps);
for sym_i = 0:N_syms-1

    datasym = sig(datasym1_start_ind + sym_i*sym_len + GI_len : datasym1_start_ind + sym_i*sym_len + GI_len +   fft_len - 1);
    datasym_fft = fft(datasym)/fft_len;
    qam_mod_syms_scaled = datasym_fft(inds)*sqrt(56)/K_mod; % scaled for qamdemod
    qam_syms = qamdemod(qam_mod_syms_scaled,M);
    a = fliplr(de2bi(qam_syms));
    bits = reshape(a.',1,numel(a));
%     bits_out = horzcat(bits_out, bits);
    bits_out(sym_i*N_dbps+1:(sym_i+1)*N_dbps) = bits;
end

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
