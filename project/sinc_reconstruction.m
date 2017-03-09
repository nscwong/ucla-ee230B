function signal_mod = sinc_reconstruction(signal, tau_t, Ts)
% Recreate the signal into signal_mod according to the delay tau_t
% 
% signal - original signal
% tau_t - delay of signal in time
% Ts - data rate

% Filter is 201 taps total in this config
numTaps = 100;
sincTaps = -numTaps:numTaps+1;

% tau_t = 0; % tau_t is in time
tau_s = tau_t/Ts; % tau_s is in samples
signal_mod = conv(signal, sinc(-tau_s + sincTaps), 'same');

end
