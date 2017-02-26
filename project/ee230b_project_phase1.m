%% EE 230B Project
%  EWC proposal to the 802.11 standards committee
%  
%  Only Green Field mode will be considered
%  The system under study will only utilize a subset of the modes 
%  described in the standards document in Table 1. 
%    The following antenna configurations: 1x1, 1x2, and 2x2
%    spatial multiplexing
%    Only QPSK, 16-QAM and 64-QAM modulations will be studied
%    Only the 20 MHz mode will be implemented
%    Always assume a 0.25 KByte packet payload
%  No coding (convolutional or LDPC) will be included in our study

%% Initialization Parameters
f_c = 2.4e9;            % Carrier Frequency for 802.11
bandwidth = 20e6;       % 
d = 140;                % Distance between TX and RX (m)
nbits = 4;              % Data length in symbols
SNR_dB = 2;             % USER-SPECIFIED SNR

% Modulation Parameters
rrc_beta = 0;
rrc_nsym = 8;
% We need 10ns increments at minimum so we can calculate this by
% 10ns = Ts/rrc_sampspersym where Ts = 1/bandwidth
rrc_sampspersym = 5; 

% Channel Parameters
alpha = 3;              % Variable for path gain

%% Simulation

% MODULATION
Payload = randi([0 1], nbits, 1);
%Payload = [1 1 1 1 1 1 0 0 0 0 0 0 0 0]';
Payload_syms = 2*Payload - 1; % For testing purposes

rctFilt = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          rrc_beta, ...
  'FilterSpanInSymbols',    rrc_nsym, ...
  'OutputSamplesPerSymbol', rrc_sampspersym);
%fvtool(rctFilt, 'Analysis', 'impulse')

TX = step(rctFilt, [Payload_syms; zeros(rrc_nsym,1)]);

% Frequency Selective Channel:  Two Ray Model
% 50 ns difference between two multipaths
%TX_c = [TX; zeros(50,1)] + [zeros(50,1); TX];

% Frequency Selective Channel:  Four Ray Model
TX_c = [TX; zeros(200,1)] + [zeros(70,1); TX; zeros(130,1)] + ...
       [zeros(150,1); TX; zeros(50,1)] + [zeros(200,1); TX];

figure;
plot(TX, 'LineWidth', 2);
hold on
plot(TX_c, '--', 'LineWidth', 2);
hold off

% CHANNEL

% Path Gain
lambda = 3e8/f_c;
Path_Gain = (lambda/(4*pi))^2 * (1/d)^alpha;

% Shadow Fading
Shadow_dB = normrnd(0, 1, numel(TX_c), 1);
Shadow = 10.^(Shadow_dB/10);

if nbits > 10e3
figure; % PDF of log-normal shadow fading
histogram(Shadow, 'Normalization', 'pdf');
title('Shadow Fading Channel Model: Amplitude');
xlabel('Amplitude of Shadow');
ylabel('Probability Density Function');

figure; % PDF of log-normal shadow fading
histogram(10*log10(Shadow), 'Normalization', 'pdf');
title('Shadow Fading Channel Model: Amplitude');
xlabel('Amplitude of Log-Normal Shadow');
ylabel('Probability Density Function');
end

% Rayleigh Fading
ComplexFading = normrnd(0, 1, numel(TX_c), 1) + 1i.*normrnd(0, 1, numel(TX_c), 1);
ComplexFading_abs = abs(ComplexFading);
ComplexFading_ang = atan2(imag(ComplexFading), real(ComplexFading));

if nbits > 10e3
figure; % PDF of fast fading amplitude
histogram(ComplexFading_abs, 'Normalization', 'pdf');
title('Fast Fading Channel Model: Amplitude');
xlabel('Amplitude of Complex Gaussian');
ylabel('Probability Density Function');

figure; % PDF of fast fading phase
histogram(ComplexFading_ang, (-pi):(pi/48):(pi), 'Normalization', 'pdf');
title('Fast Fading Channel Model: Phase');
xlabel('Phase of Complex Gaussian (radians)');
ylabel('Probability Density Function');
end

RX_1 = TX_c;
RX_2 = RX_1*sqrt(Path_Gain);
RX_3 = RX_2.*sqrt(Shadow);
RX_4 = RX_3.*sqrt(abs(ComplexFading));

% AWGN Noise
Signal = RX_4;
SNR = 10^(SNR_dB/10);                       % Convert to linear
Esym = sum(abs(Signal).^2)/numel(Signal);   % Actual symbol energy
N0 = Esym/SNR;
if (isreal(Signal))
    noise_sigma = sqrt(N0);
    Noise = noise_sigma*randn(numel(Signal),1);
else
    noise_sigma = sqrt(N0/2);
    Noise = noise_sigma*(randn(numel(Signal),1) + 1j*randn(numel(Signal),1));
end

RX_5 = RX_4 + Noise;

if nbits <= 50
figure;
stem(RX_2);
hold on
stem(RX_3);
stem(RX_4);
stem(RX_5);
hold off
legend('Path Gain', '+ Shadow Fading', '+ Fast Fading', '+ Noise');
end