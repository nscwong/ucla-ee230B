
string1 = 'awgn_channel';
string2 = 'rayleigh_channel';
string3 = '_ideal';
string4 = '_nonideal';
string5 = '_agc';
string6 = '_bbd';
string7 = '_chanest';
string8 = '_freqest';
string9 = '.mat';
string10 = '_extra';
string11 = '../rayleigh_m04_nomulti/';
string12 = '../rayleigh_m16_nomulti/';

figure;
hold all
string = [string12 string2 string3 string9];
load(string);
semilogy(SNR_Values, SNR_ErrorRates','b-','LineWidth',1.5);
string = [string12 string2 string4 string5 string6 string7 string8 string9];
load(string);
semilogy(SNR_Values, SNR_ErrorRates','b--','LineWidth',1.5);
string = [string2 string3 string9];
load(string);
semilogy(SNR_Values, SNR_ErrorRates','r-','LineWidth',1.5);
string = [string2 string4 string5 string6 string7 string8 string9];
load(string);
semilogy(SNR_Values, SNR_ErrorRates','r--','LineWidth',1.5);

title('BER for Rayleigh Channel');
xlabel('SNR Values');
ylabel('Bit Error Rate (BER)');
legend('Ideal Simulation -- 1x1','NonIdeal All -- 1x1','Ideal Simulation -- 1x2','NonIdeal All -- 1x2','Location','best');
