% This is a Rayleigh Distribution (We'll probably need it later)

t = 1000; % Random value?

C = normrnd(0, 1, t) + 1i.*normrnd(0, 1, t); % Is this log-normal?
C_abs = abs(C);
C_ang = atan2(imag(C), real(C));

figure; % PDF of amplitude
histogram(C_abs, 'Normalization', 'pdf');
title('Flat Fading Channel Model: Amplitude');
xlabel('Amplitude of Complex Gaussian');
ylabel('Probability Density Function');

figure; % PDF of phase
histogram(C_ang, (-pi):(pi/48):(pi), 'Normalization', 'pdf');
title('Flat Fading Channel Model: Phase');
xlabel('Phase of Complex Gaussian (radians)');
ylabel('Probability Density Function');