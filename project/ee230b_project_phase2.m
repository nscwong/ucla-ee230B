d = 100;
SNR_dB = 20.0;
amplitudes = [0, 0];
delays = [0, 50e-9];

iss = 1;
nss = 1;
tx_packet = create_tx_packet([],iss,nss);
packet = create_channel_model(tx_packet,amplitudes,delays,d,SNR_dB);
t_start = 1;
t_end = 500;
% Insert AGC here
threshold = 0.5;
[packet_detected, decision_statistics] = packet_detection(packet, iss, nss, t_start, t_end, threshold);
[test, test2] = htltf_start(packet, iss, nss, t_start, t_end);

