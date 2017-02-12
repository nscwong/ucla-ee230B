%% Reference code for RRC TX and RX code
%https://www.mathworks.com/help/comm/examples/raised-cosine-filtering.html?refresh=true#zmw57dd0e6093
Nsym = 6;               % Filter span in symbol durations
rrc_beta = 0.5;         % Roll-off factor
sampsPerSym = 8;        % Upsampling factor

% Parameters
DataL = 20;             % Data length in symbols
R = 40e6;               % Data rate
Fs = R * sampsPerSym;   % Sampling frequency

% Create a local random stream to be used by random number generators for
% repeatability
hStr = RandStream('mt19937ar', 'Seed', 0);

% Generate random data
x = 2*randi(hStr, [0 1], DataL, 1)-1;

% Time vector sampled at symbol rate in milliseconds
time = 1000 * (0: DataL - 1) / R;

% Design raised cosine filter with given order in symbols
rctFilt3 = comm.RaisedCosineTransmitFilter( ...
    'Shape',                    'Square root', ...
    'RolloffFactor',            rrc_beta, ...
    'FilterSpanInSymbols',      Nsym, ...
    'OutputSamplesPerSymbol',   sampsPerSym);

%fvtool(rctFilt3, 'Analysis', 'impulse')

% Upsample and filter.
yc = rctFilt3([x; zeros(Nsym/2,1)]);

% Filter group delay, since raised cosine filter is linear phase and
% symmetric.
fltDelay = Nsym / (2*R);
% Correct for propagation delay by removing filter transients
yc = yc(fltDelay*Fs+1:end);

% Plot data.
figure(1)
stem(time, x, 'kx'); hold on;
% Plot filtered data.
to = 1000 * (0: DataL*sampsPerSym - 1) / Fs;
plot(to, yc, 'm-'); hold off;
% Set axes and labels.
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Transmitted Data', 'Sqrt. Raised Cosine', 'Location', 'southeast')


% Design and normalize filter.
rcrFilt = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          rrc_beta, ...
  'FilterSpanInSymbols',    Nsym, ...
  'InputSamplesPerSymbol',  sampsPerSym, ...
  'DecimationFactor',       1);
% Filter at the receiver.
yr = rcrFilt([yc; zeros(Nsym*sampsPerSym/2, 1)]);
% Correct for propagation delay by removing filter transients
yr = yr(fltDelay*Fs+1:end);
% Plot data.
figure(2)
stem(time, x, 'kx'); hold on;
% Plot filtered data.
plot(to, yr, 'b-'); hold off;
% Set axes and labels.
xlabel('Time (ms)'); ylabel('Amplitude');
legend('Transmitted Data', 'Rcv Filter Output',...
    'Location', 'southeast')