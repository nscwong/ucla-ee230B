function tx_packet = create_tx_packet(data)
% Creates a TX Packet Structure according to the EWC standard
% 
% Ignore: Encoder Parser & FEC Encoder
% Remember: Only the 20 MHz mode will be implemented

% Create the subframes
% Page 
subframe_lstf = sqrt(0.5)*[0,0,1+1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,0,0,0,0,-1-1j,0,0,0,-1-1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0,0,1+1j,0,0];
% Page 26 and 27 of spec
subframe_htltf1 = [1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,1,1,-1,-1,1,1,-1,1,-1,1,1,1,1,0,1,-1,-1,1,1,-1,1,-1,1,-1,-1,-1,-1,-1,1,1,-1,-1,1,-1,1,-1,1,1,1,1,-1,-1];
% Page 21 and 22 of spec
mcs = zeros(1,7);       % Modulation and Coding Scheme Table 1 LSB First
bw2040 = 0;             % 0 if 20MHz or 40MHz upper/lower, 1 if 40MHz
length = zeros(1,16);   % The number of bytes in frame - LSB first
aggregation = 0;        % 1 if PPDU in data packet contains an A-MPDU
stbc = zeros(1,2);      %
advanced_coding = 0;    % 1 if Advanced Coding, 0 if BCC
short_gi = 0;           % Indicate that short GI is used after HT training
nhtltf = zeros(1,2);    % Number of HT-LTF
crc = zeros(1,7);       % CRC bits 
tail_bits = zeros(1,6); % Tail bits
subframe_htsig = [mcs,bw2040,length,1,1,1,aggregation,stbc,advanced_coding,short_gi,nhtltf,crc,tail_bits];
% Page 27 and 28 of spec
P_htltf = [ ...     % Polarity matrix
     1 -1  1  1;    % P_htlf(i,n) represents polarity of the
     1  1 -1  1;    % ith spatial stream in the
     1  1  1 -1;    % nth HT training symbol
    -1  1  1  1;
];
subframe_htltfs = [];
% Page 28 and 29 of spec
subframe_data = data;

tx_packet = [subframe_lstf, subframe_htltf1, subframe_htsig, subframe_htltfs, subframe_data];

% Plot data


end