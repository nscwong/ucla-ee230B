function tx_packets = create_tx_packet(data, nss)
% Creates a TX Packet Structure according to the EWC standard
% 
% nss - Number of Spatial Streams
% iss - This spatial stream ID
%
% Ignore: Encoder Parser & FEC Encoder
% Remember: Only the 20 MHz mode will be implemented

% To replace create_tx_packet

tx_packets = cell(nss,1);

% Initial parameters -- copied for convenience
f_c = 2.4e9;            % Carrier Frequency for 802.11
bandwidth = 20e6;       % Bandwidth
dt = 1/bandwidth;       % Sample Period
df = 312.5e3;           % Subcarrier frequency spacing (delta f)
nbits = 4;              % Data length in symbols

M = 4;      % M-QAM number of modulation points
data = randi([0 1], 1, 9);  % 1xn vector of data

nfft = 64;

% Timing related constants -- Page 11 of spec, Table 2
NSD_l = 48;                 % Number of data subcarriers for legacy
NSP_l = 4;                  % Number of pilot subcarriers for legacy
NST_l = NSD_l + NSP_l;      % Total number of subcarriers for legacy
NSR_l = 26;                 % Num subcarriers occupying half of overall BW in legacy
NSD_ht = 52;                % NSD for HT
NSP_ht = 4;                 % NSP for HT
NST_ht = NSD_ht + NSP_ht;   % NST for HT
NSR_ht = 28;                % NSR for HT
T_lstf = 8e-6;
T_htltf1 = 8e-6;
T_htltfs = 4e-6;
T_GI = 0.8e-6;

time = 0:dt:(T_lstf+T_htltf1+3*T_htltfs)-dt; % Time vector for Preamble

% Cyclic shift for Legacy Preamble
T_iTX_CS_l = [ ...
    0    0    0    0;
    0 -200    0    0;
    0 -100 -200    0;
    0  -50 -100 -150;
];
% Cyclic shift for High Throughput (HT) Preamble -- Page 20 of spec
T_iTX_CS_ht = [ ...
    0    0    0    0;
    0 -400    0    0;
    0 -400 -200    0;
    0 -400 -200 -600;
];
% Polarity matrix -- Page 27 and 28 of spec
P_htltf = [ ...     % 
     1 -1  1  1;    % P_htlf(i,n) represents polarity of the
     1  1 -1  1;    % ith spatial stream in the
     1  1  1 -1;    % nth HT training symbol
    -1  1  1  1;    %
];

% Number of tones in each field -- Page 15 of spec, Table 4
Ntone_lstf = 12;
Ntone_htsig = 56;
Ntone_htltf = 56;
Ntone_htdata = 56;

for iss = 1:nss

    % L-STF -- Page
    symmap_lstf = sqrt(0.5)*[0,0,1+1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,0,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0];
    t_lstf = 0:dt:(T_lstf-dt);
    subframe_lstf = zeros(size(t_lstf));
    for k = -NSR_l:NSR_l
        subframe_lstf = subframe_lstf + symmap_lstf(k+NSR_l+1)*...
                        exp(2j*pi*k*df*(t_lstf-T_iTX_CS_l(nss,iss)));
    end
    subframe_lstf = subframe_lstf/sqrt(Ntone_lstf);

    % HT-LTF1 -- Page 26 to 28 of spec
    symmap_htltf = [1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1];
    t_htltf1 = 0:dt:(T_htltf1-dt);
    subframe_htltf1 = zeros(size(t_htltf1));
    for k = -NSR_ht:NSR_ht
        subframe_htltf1 = subframe_htltf1 + ...
                          P_htltf(iss,1)*symmap_htltf(k+NSR_ht+1)*...
                          exp(2j*pi*k*df*(t_htltf1-2*T_GI-T_iTX_CS_ht(nss,iss)));
    end
    subframe_htltf1 = subframe_htltf1/sqrt(Ntone_htltf);

    % NOTE: We omit HT-SIG for this project
    % Page 21 and 22 of spec
    % mcs = zeros(1,7);       % Modulation and Coding Scheme Table 1 LSB First
    % bw2040 = 0;             % 0 if 20MHz or 40MHz upper/lower, 1 if 40MHz
    % length = zeros(1,16);   % The number of bytes in frame - LSB first
    % aggregation = 0;        % 1 if PPDU in data packet contains an A-MPDU
    % stbc = zeros(1,2);      %
    % advanced_coding = 0;    % 1 if Advanced Coding, 0 if BCC
    % short_gi = 0;           % Indicate that short GI is used after HT training
    % nhtltf = zeros(1,2);    % Number of HT-LTF
    % crc = zeros(1,7);       % CRC bits 
    % tail_bits = zeros(1,6); % Tail bits
    % subframe_htsig = [mcs,bw2040,length,1,1,1,aggregation,stbc,advanced_coding,short_gi,nhtltf,crc,tail_bits];

    % HT-LTFs -- Page 27 to 28 of spec
    t_htltfs = 0:dt:(T_htltfs-dt);
    subframe_htltfs = zeros(size(t_htltfs));
    for k = -NSR_ht:NSR_ht
        subframe_htltfs = subframe_htltfs + symmap_htltf(k+NSR_ht+1)*...
                          exp(2j*pi*k*df*(t_htltfs-2*T_GI-T_iTX_CS_ht(nss,iss)));
    end
    subframe_htltfs = subframe_htltfs/sqrt(Ntone_htltf);

    subframe_htltf2 = P_htltf(iss,2)*subframe_htltfs;
    subframe_htltf3 = P_htltf(iss,3)*subframe_htltfs;
    subframe_htltf4 = P_htltf(iss,4)*subframe_htltfs;
    
    tx_packets{iss} = [subframe_lstf,subframe_htltf1,subframe_htltf2,subframe_htltf3,subframe_htltf4];
end
% Data
padamt = M - mod(numel(data),M);    % Number of zeros to pad
data = [data zeros(1,padamt)];      % Pad to fit modulation scheme
data_i = reshape(data,log2(M),numel(data)/log2(M));    % Change to int
data_i = (2.^((log2(M)-1):-1:0))*data_i;
data_mod = qammod(data_i,M);  % Modulated data

% FIXME:  Fix Data Packet

% Final Packet

% Plot first signal
figure;
plot(time, real(tx_packets{1}));
title('create\_tx\_packet');

end