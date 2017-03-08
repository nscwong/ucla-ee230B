function [packet_detected, decision_statistics] = packet_detection(tx_packet, iss, nss, t_start, t_end, threshold)
% Packet Detection Function
% 
% Outputs:
% Returns packet_detected as 1 if packet detected, 0 if not
% Returns decision_statistics the peak convolution
%
% Inputs:
% tx_packet - transmitted packet
% iss - stream ID
% nss - number of spatial streams
% t_start - when to start looking at packet
% t_end - when to stop looking at packet
% threshold - when a peak is determined in decision_statistics

% We assume that 2 peaks are sufficient
PEAKS_NEEDED = 2;
% We assume that the distance between the first and subsequent peaks
% must be less than 3 Ts
PERIOD_LIMIT = 3;

packet_detected = 0;

% Initial parameters -- copied for convenience
f_c = 2.4e9;            % Carrier Frequency for 802.11
bandwidth = 20e6;       % Bandwidth
df = 312.5e3;           % Subcarrier frequency spacing (delta f)

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

% Number of tones in each field -- Page 15 of spec, Table 4
Ntone_lstf = 12;

% Create 1 period worth of the STF signal
symmap_lstf = sqrt(0.5)*[0,0,1+1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,0,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0];
dt = 1/bandwidth;
t_lstf = 0:dt:(T_lstf-dt);
t_lstf_compare = 0:dt:0.8e-6-dt;
compare_lstf = zeros(size(t_lstf));
for k = -NSR_l:NSR_l
    compare_lstf = compare_lstf + symmap_lstf(k+NSR_l+1)*...
                    exp(2j*pi*k*df*(t_lstf-T_iTX_CS_l(nss,iss)));
end
compare_lstf = compare_lstf/sqrt(Ntone_lstf);
compare_lstf = compare_lstf(1:numel(t_lstf_compare));
% figure;
% plot(t_lstf_compare, compare_lstf);

% Option 1: Correlation Attempt
% c_1 = xcorr(tx_packet(t_start:t_end),compare_lstf);
% p_1 = sum(abs(compare_lstf).^2);
% m_1 = abs(c_1).^2/p_1^2;
% figure;
% plot(m_1);

% Option 2: Convolution Attempt
c_2 = conv(tx_packet(t_start:t_end), compare_lstf);
p_2 = sum(abs(compare_lstf).^2);
m_2 = abs(c_2).^2/p_2^2;
figure;
plot(m_2(t_start:t_end));

% Option 3: Periodicity Attempt
% https://www.mathworks.com/help/wlan/ref/wlanpacketdetect.html#bvaykqk
% period = 0:dt:0.8e-6;
% if t_start > numel(period)
%     c_3 = xcorr(tx_packet(t_start:t_end-1), tx_packet(t_start-numel(period):t_start-1));
%     p_3 = sum(abs(tx_packet(t_start-numel(period):t_start-1)).^2);
%     m_3 = abs(c_3).^2/p_3^2;
%     figure;
%     plot(m_3);
% end

m = m_2;
% Count the number of peaks and the respective distances
npeaks = sum(m > threshold);
distance = find(m > threshold, 2);
if (npeaks > PEAKS_NEEDED) & (distance < PERIOD_LIMIT*numel(t_lstf))
    packet_detected = 1;
end
decision_statistics = m(t_start:t_end);
