function tx_packet = create_tx_packet()
% Creates a TX Packet Structure according to the EWC standard
% 
% Ignore: Encoder Parser & FEC Encoder
% Remember: Only the 20 MHz mode will be implemented

% Create the subframes
subframe_lstf = sqrt(0.5)*[0,0,1+1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,0,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0];
subframe_htltf1 = [];
subframe_htsig = [];
subframe_htltfs = [];
subframe_data = [];

end