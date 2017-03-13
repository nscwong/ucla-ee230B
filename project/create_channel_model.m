function packet = create_channel_model(tx_packet, amplitudes_dB, delays, d, SNR_dB)
% Channel Model for one packet
%
% d - Distance between TX and RX (m)
% SNR_dB - USER-SPECIFIED SNR

% Initialization Parameters
f_c = 2.4e9;            % Carrier Frequency for 802.11
bandwidth = 20e6;       %  
dt = 1/bandwidth;

nmultipath = numel(amplitudes_dB);
amplitudes = 10.^(amplitudes_dB/10); % Convert from dB

% Create Path Loss/Gain parameter
alpha = 3;          % Path loss exponent
lambda = 3e8/f_c;   % wavelength, c = 3e8 m/s speed of light
Path_Gain = (lambda/(4*pi))^2 * (1/d)^alpha;

% Create Shadowing parameter - one value per packet
Shadow_dB = normrnd(0, 1);
Shadow = 10.^(Shadow_dB/10);

% Create Rayleigh Fading parameter - one value per multipath component
Complex_Fading = normrnd(0, 1, nmultipath, 1) ...
               + 1i.*normrnd(0, 1, nmultipath, 1);

% Create Noise component
Signal = zeros(size(tx_packet));
for k = 1:nmultipath
    mp = create_multipath(tx_packet, amplitudes(k), delays(k), dt);
    Signal = Signal + mp;
end
Signal = Signal*Path_Gain;
SNR = 10^(SNR_dB/10);                       % Convert to linear
Esym = sum(abs(Signal).^2)/numel(Signal);   % Actual symbol energy
N0 = Esym/SNR;
if (isreal(Signal))
    noise_sigma = sqrt(N0);
    Noise = noise_sigma*randn(1,numel(Signal));
else
    noise_sigma = sqrt(N0/2);
    Noise = noise_sigma*(randn(1,numel(Signal)) + 1j*randn(1,numel(Signal)));
end

sig_0 = zeros(size(tx_packet)); % Only multipath
sig_1 = sig_0;                  % Add path gain
sig_2 = sig_0;                  % Add shadowing
sig_3 = sig_0;                  % Add complex fading
sig_4 = sig_0;                  % Add noise

for k = 1:nmultipath
    mp = create_multipath(tx_packet, amplitudes(k), delays(k), dt);
    sig_0 = sig_0 + mp;
    mp = mp*Path_Gain;
    sig_1 = sig_1 + mp;
    mp = mp*Shadow;
    sig_2 = sig_2 + mp;
    mp = mp*Complex_Fading(k);
    sig_3 = sig_3 + mp;
    mp = mp + Noise;
    sig_4 = sig_4 + mp;
end

packet = sig_4;

time_x = 0:dt:(numel(tx_packet)*dt-dt);

figure;
plot(time_x, real(sig_1));
hold on
plot(time_x, real(sig_2));
plot(time_x, real(sig_3));
plot(time_x, real(sig_4));
hold off
title('Packet in Channel (Real)');
legend('Path Gain', '+ Shadowing', '+ Complex Fading', '+ Noise');

figure;
plot(time_x, imag(sig_1));
hold on
plot(time_x, imag(sig_2));
plot(time_x, imag(sig_3));
plot(time_x, imag(sig_4));
hold off
title('Packet in Channel (Imaginary)');
legend('Path Gain', '+ Shadowing', '+ Complex Fading', '+ Noise');

end