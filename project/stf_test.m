% Initial parameters -- copied for convenience
f_c = 2.4e9;            % Carrier Frequency for 802.11
bandwidth = 20e6;       % Bandwidth
df = 312.5e3;           % Subcarrier frequency spacing (delta f)
nbits = 4;              % Data length in symbols

nss = 1;    % Number of spatial streams
iss = 1;    % This spatial stream ID

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

% L-STF -- Page
symmap_lstf = sqrt(0.5)*[0,0,1+1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,0,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0];
t_lstf = 0:(1/bandwidth):T_lstf;
subframe_lstf = zeros(size(t_lstf));
for k = -NSR_l:NSR_l
    subframe_lstf = subframe_lstf + symmap_lstf(k+NSR_l+1)*...
                    exp(2j*pi*k*df*(t_lstf-T_iTX_CS_l(nss,iss)));
end
subframe_lstf = subframe_lstf/sqrt(Ntone_lstf);

% HT-LTF1 -- Page 26 to 28 of spec
symmap_htltf = [1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1];
t_htltf1 = 0:(1/bandwidth):T_htltf1;
subframe_htltf1 = zeros(size(t_htltf1));
for k = -NSR_ht:NSR_ht
    subframe_htltf1 = subframe_htltf1 + ...
                      P_htltf(iss,1)*symmap_htltf(k+NSR_ht+1)*...
                      exp(2j*pi*k*df*(t_lstf-2*T_GI-T_iTX_CS_l(nss,iss)));
end
subframe_htltf1 = subframe_htltf1/sqrt(Ntone_htltf);

% Final Packet
time = [t_lstf,(t_htltf1+T_lstf)];
signal = [subframe_lstf,subframe_htltf1];

% Plot signal
figure;
plot(time, real(signal));
