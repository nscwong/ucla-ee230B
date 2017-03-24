function packet = create_channel_model(tx_packet, powers_dB, delays, d, SNR_dB)
% Channel Model for one packet
%
% d - Distance between TX and RX (m)
% SNR_dB - USER-SPECIFIED SNR

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

Complex_Fading = ones(nmultipath,numel(tx_packet));
% for k = 1:nmultipath
%     % Create Rayleigh Fading parameter - one set per multipath component
%     C = normrnd(0, 1, numel(tx_packet), 1) + 1i.*normrnd(0, 1, numel(tx_packet), 1);
%     Fs = bandwidth;
%     half_shift = floor(length(C)/2);
%     C_fft_shift = circshift(fft(C)/length(C),-half_shift);
%     % fft frequency list (X-axis)
%     Freqs = 0:(length(C)-1);
%     Freqs = Freqs - half_shift;
%     Freqs = Freqs*Fs/length(C);
%     Fdoppler = 200; % Maximum doppler frequency
%     % Calculate PSD
%     PSD = 2/(pi*Fdoppler)./sqrt(1-(Freqs/Fdoppler).^2);
%     PSD(abs(Freqs)>=Fdoppler) = 0; % Zero out PSD outside of valid range
%     C_fft_shift_filtered = C_fft_shift.*sqrt(PSD'); %Apply filter
%     C_filtered = ifft(circshift(C_fft_shift_filtered,half_shift));
%     Complex_Fading(k,:) = C_filtered';
%     figure;
%     plot(real(Complex_Fading(k,:)));
%     hold on
%     plot(imag(Complex_Fading(k,:)));
%     hold off
% end

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
        figure;
        plot(Rayleigh_Channel{k});
    else
        mp = mp*Complex_Fading(k);
    end
    sig_3 = sig_3 + mp;    
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

ExtraNoise = noise_sigma*(randn(nmultipath,300) + 1j*randn(nmultipath,300));
if nmultipath > 1
    ExtraNoise = sum(ExtraNoise);
end

packet = sig_4;
packet = [ExtraNoise sig_4];

time_x = 0:dt:(numel(tx_packet)*dt-dt);

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