function tx_packet_multipath = create_multipath(tx_packet, amplitude, delay, Ts)
% Create a single multipath component for tx_packet

tx_packet_multipath = amplitude*sinc_reconstruction(tx_packet, delay, Ts);
end
