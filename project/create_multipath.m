function tx_packet_multipath = create_multipath(tx_packet, power, delay, Ts)
% Create a single multipath component for tx_packet

tx_packet_multipath = sqrt(power)*sinc_reconstruction(tx_packet, delay, Ts);
end
