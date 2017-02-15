d = [60 100 140];
t_outage = [0.12565  0.6973   0.94659];
s_outage = [0.002083 0.014813 0.047417];

figure;
plot(d, t_outage, '-o');
hold on
plot(d, s_outage, '-o');
hold off
title('Outage of System');
xlabel('Distance between TX and RX (d)');
ylabel('Outage');
legend('Theoretical', 'Simulation', 'Location', 'best');