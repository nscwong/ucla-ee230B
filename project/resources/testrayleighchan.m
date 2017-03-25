bandwidth = 1e6;
h = rayleighchan(1/bandwidth, 10);
hMod = comm.DPSKModulator('ModulationOrder',2);
tx = randi([0 1],1e6,1);          % Random bit stream
dpskSig = step(hMod,tx);          % DPSK signal
h.StoreHistory = true;            % Allow states to be stored
y = filter(h, dpskSig);           % Run signal through channel
plot(h);                          % Call Channel Visualization Tool