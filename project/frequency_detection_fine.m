function fine_freq_offset_est_hz = frequency_detection_fine(packet_fc_coarse, start_ind_htltf, coarse_freq_offset_est_hz)

Fs = 20e6; % Sample Frequency
Fc = 2.4e9; % Carrier Frequency
LTF_len = 64;
LTF_sig = packet_fc_coarse(start_ind_htltf:start_ind_htltf+LTF_len*2-1);

delta_phase = phase(LTF_sig(LTF_len+1:LTF_len*2)*LTF_sig(1:LTF_len)');
fine_freq_offset_est_hz = delta_phase/(2*pi*LTF_len/Fs);
% fine_freq_offset_est_ppm = fine_freq_offset_est_hz/Fc*1e6;
% 
% disp(['Fine Freq Offset Est (ppm): ', num2str(fine_freq_offset_est_ppm)]);

% total_freq_offset_est_hz = coarse_freq_offset_est_hz+fine_freq_offset_est_hz;
% total_freq_offset_est_ppm = total_freq_offset_est_hz/Fc*1e6;
% 
% disp(['Total Freq Offset Est (ppm): ', num2str(total_freq_offset_est_ppm)]);

end