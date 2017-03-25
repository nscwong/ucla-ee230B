% IN PROGRESS
function [start_point, decision_statistics] = htltf_start(tx_packet, iss, nss, t_start, t_end, OutputPlots)

if nargin < 6
    OutputPlots = 0;
end

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

% Create 1 period worth of the STF signal
symmap_lstf = sqrt(0.5)*[0,0,1+1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,0,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0];
dt = 1/bandwidth;
t_lstf = 0:dt:(T_lstf-dt);
t_lstf_compare = 0:dt:4e-6-dt;
%t_lstf_compare = 0;
compare_lstf = zeros(size(t_lstf));
for k = -NSR_l:NSR_l
    compare_lstf = compare_lstf + symmap_lstf(k+NSR_l+1)*...
                   exp(2j*pi*k*df*(t_lstf-T_iTX_CS_l(nss,iss)));
end
compare_lstf = compare_lstf/sqrt(Ntone_lstf);
compare_lstf = compare_lstf(end-numel(t_lstf_compare)+1:end);
%compare_lstf = [];

% Create 1 period worth of the HT-LTF1
symmap_htltf = [1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1];
t_htltf1 = 0:dt:(T_htltf1-dt);
t_htltf1_compare = 0:dt:(4e-6-dt);
compare_htltf1 = zeros(size(t_htltf1));
for k = -NSR_ht:NSR_ht
    compare_htltf1 = compare_htltf1 + ...
                     P_htltf(iss,1)*symmap_htltf(k+NSR_ht+1)*...
                     exp(2j*pi*k*df*(t_htltf1-2*T_GI-T_iTX_CS_ht(nss,iss)));
end
compare_htltf1 = compare_htltf1/sqrt(Ntone_htltf);
compare_htltf1 = compare_htltf1(1:numel(t_htltf1_compare));

% c_l = filter(fliplr(conj(compare_lstf)),1,tx_packet(t_start:t_end));
% p_l = sum(abs(compare_lstf).^2);
% m_l = abs(c_l).^2/p_l^2;
% figure;
% plot(m_l(numel(t_lstf_compare):end));
% 
% c_h = filter(fliplr(conj(compare_htltf1)),1,tx_packet(t_start:t_end));
% p_h = sum(abs(compare_htltf1).^2);
% m_h = abs(c_h).^2/p_h^2;
% figure;
% plot(m_h(numel(t_htltf1_compare)+numel(t_lstf):end));

%compare_filter = fliplr(conj([compare_lstf compare_htltf1]));
compare_filter = fliplr(conj([compare_lstf compare_htltf1]));
c = filter(compare_filter,1,tx_packet(t_start:t_end));
p = sum(abs(compare_filter).^2);
% p = sum(abs(tx_packet(t_start:t_end)).^2);
% p(isnan(p)) = 1;
m = abs(c).^2/p^2;

if OutputPlots
figure;
plot(m)
title('HT-LTF Start');
end

% figure;
% plot(real(compare_filter));
% hold on
% plot(imag(compare_filter));
% hold off

decision_statistics = m;
[~,start_point] = max(m);
start_point = start_point - numel(t_htltf1_compare) + t_start - 1;

end