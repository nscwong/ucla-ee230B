% This is a Rayleigh Distribution (We'll probably need it later)

t = 1e5; % Random value?

C = normrnd(0, 1, t, 1) + 1i.*normrnd(0, 1, t, 1); 
C_abs = abs(C);
C_ang = atan2(imag(C), real(C));

% figure; % PDF of amplitude
% histogram(C_abs, 'Normalization', 'pdf');
% title('Flat Fading Channel Model: Amplitude');
% xlabel('Amplitude of Complex Gaussian');
% ylabel('Probability Density Function');
% 
% figure; % PDF of phase
% histogram(C_ang, (-pi):(pi/48):(pi), 'Normalization', 'pdf');
% title('Flat Fading Channel Model: Phase');
% xlabel('Phase of Complex Gaussian (radians)');
% ylabel('Probability Density Function');

Fs = 20e6;
half_shift = floor(length(C)/2);
C_fft_shift = circshift(fft(C)/length(C),-half_shift);

figure(1)
% fft frequency list (X-axis)
Freqs = 0:(length(C)-1);
Freqs = Freqs - half_shift;
Freqs = Freqs*Fs/length(C);
plot(Freqs,10*log10(abs(C_fft_shift)),'.')
% hold all

% figure(2)
% pwelch(C,[],[],[],Fs);

Fdoppler = 100; % Maximum doppler frequency

% Calculate PSD
PSD = 2/(pi*Fdoppler)./sqrt(1-(Freqs/Fdoppler).^2);
PSD(abs(Freqs)>=Fdoppler) = 0; % Zero out PSD outside of valid range
C_fft_shift_filtered = C_fft_shift.*sqrt(PSD'); %Apply filter

figure(2)
% plot(PSD)

plot(Freqs,10*log10(abs(C_fft_shift_filtered)),'.-')
plot(Freqs,abs(C_fft_shift_filtered),'.-')
xlim([-300,300])


% Undue the circular shift and ifft to get time domain signal
C_filtered = ifft(circshift(C_fft_shift_filtered,half_shift));

time = 1:t;
time = time./Fs;

figure(3)
plot(time,C_filtered,'.-')

figure(4)
pwelch(C_filtered,[],[],[],Fs)

