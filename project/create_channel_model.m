function [packet, packet_no_awgn] = create_channel_model(tx_packet, powers_dB, delays, d, SNR_dB, AWGN_Only, OutputPlots)
% Channel Model for one packet
%
% d - Distance between TX and RX (m)
% SNR_dB - USER-SPECIFIED SNR

if nargin < 6
    AWGN_Only = 0;
end
if nargin < 7
    OutputPlots = 0;
end

RAYLEIGH_MATLAB = 1;

% Initialization Parameters
f_c = 2.4e9;            % Carrier Frequency for 802.11
bandwidth = 20e6;       %  
dt = 1/bandwidth;

nmultipath = numel(powers_dB);
powers = 10.^(powers_dB/10); % Convert from dB

% Create Path Loss/Gain parameter
alpha = 3;          % Path loss exponent
lambda = 3e8/f_c;   % wavelength, c = 3e8 m/s speed of light
Path_Gain = (lambda/(4*pi))^2 * (1/d)^alpha;

% Create Shadowing parameter - one value per packet
shadow_stddev_dB = 4;
Shadow_dB = normrnd(0, shadow_stddev_dB, 1, nmultipath);
Shadow = 10.^(Shadow_dB/10);

% Create Rayleigh Fading parameter - one set per multipath component
Complex_Fading = ones(nmultipath,numel(tx_packet));
Fdoppler = 10; % Maximum doppler frequency
if RAYLEIGH_MATLAB
    Rayleigh_Channel = cell(1,nmultipath);
    for k = 1:nmultipath
        Rayleigh_Channel{k} = rayleighchan(1/bandwidth,Fdoppler);
        Rayleigh_Channel{k}.StoreHistory = true;
    end
end

sig_0 = zeros(size(tx_packet)); % Only multipath
sig_1 = sig_0;                  % Add path gain
sig_2 = sig_0;                  % Add shadowing
sig_3 = sig_0;                  % Add complex fading               

for k = 1:nmultipath
    mp = create_multipath(tx_packet, powers(k), delays(k), dt);
    sig_0 = sig_0 + mp;
    mp = mp*Path_Gain;
    sig_1 = sig_1 + mp;
    mp = mp*Shadow(k);
    sig_2 = sig_2 + mp;
    if RAYLEIGH_MATLAB
        mp = filter(Rayleigh_Channel{k},mp);
    else
        mp = mp*Complex_Fading(k);
    end
    sig_3 = sig_3 + mp;    
end
if AWGN_Only
   sig_3 = sig_1; 
end

% Create Noise component
Signal = zeros(size(tx_packet));
for k = 1:nmultipath
    mp = create_multipath(tx_packet, powers(k), delays(k), dt);
    mp = mp*Path_Gain;
    if ~RAYLEIGH_MATLAB
        mp = mp*Complex_Fading(k);
    end
    Signal = Signal + mp;
end
SNR = 10^(SNR_dB/10);                       % Convert to linear
Esym = sum(abs(Signal).^2)/numel(Signal);   % Actual symbol energy
N0 = Esym/SNR;
if (isreal(Signal))
    noise_sigma = sqrt(N0);
    Noise = noise_sigma*randn(nmultipath,numel(Signal));
else
    noise_sigma = sqrt(N0/2);
    Noise = noise_sigma*(randn(nmultipath,numel(Signal)) + 1j*randn(nmultipath,numel(Signal)));
end
if nmultipath > 1
    Noise = sum(Noise);
end

sig_4 = sig_3 + Noise;  % Add noise

ExtraNoise1 = noise_sigma*(randn(nmultipath,300) + 1j*randn(nmultipath,300));
if nmultipath > 1
    ExtraNoise1 = sum(ExtraNoise1);
end
ExtraNoise2 = noise_sigma*(randn(nmultipath,300) + 1j*randn(nmultipath,300));
if nmultipath > 1
    ExtraNoise2 = sum(ExtraNoise2);
end

% packet = sig_4;
packet = [ExtraNoise1 sig_4 ExtraNoise2];
packet_no_awgn = [zeros(size(ExtraNoise1)) sig_3 zeros(size(ExtraNoise2))];

time_x = 0:dt:(numel(tx_packet)*dt-dt);

if OutputPlots
% figure;
% plot(time_x, real(sig_1));
% hold on
% plot(time_x, real(sig_2));
% plot(time_x, real(sig_3));
% plot(time_x, real(sig_4));
% hold off
% title('Packet in Channel (Real)');
% legend('Path Gain', '+ Shadowing', '+ Complex Fading', '+ Noise');
% 
% figure;
% plot(time_x, imag(sig_1));
% hold on
% plot(time_x, imag(sig_2));
% plot(time_x, imag(sig_3));
% plot(time_x, imag(sig_4));
% hold off
% title('Packet in Channel (Imaginary)');
% legend('Path Gain', '+ Shadowing', '+ Complex Fading', '+ Noise');

figure;
plot(time_x, abs(sig_1));
hold on
plot(time_x, abs(sig_2));
plot(time_x, abs(sig_3));
plot(time_x, abs(sig_4));
hold off
title('Packet in Channel (Abs)');
legend('Path Gain', '+ Shadowing', '+ Complex Fading', '+ Noise');

figure;
plot(time_x, phase(sig_1));
hold on
plot(time_x, phase(sig_2));
plot(time_x, phase(sig_3));
plot(time_x, phase(sig_4));
hold off
title('Packet in Channel (Abs)');
legend('Path Gain', '+ Shadowing', '+ Complex Fading', '+ Noise');

% time_xmore = 0:dt:(numel(packet)*dt-dt);
% figure;
% plot(real(packet));
% hold on
% plot(imag(packet));
% hold off
end

end