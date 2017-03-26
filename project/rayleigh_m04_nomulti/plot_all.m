
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

figure;
hold all
for i = 1:5
    if i == 1
        string = [string2 string3 string9];
    end
    if i == 2
        string = [string2 string4 string5 string9];
        string_extra = [string1 string4 string5 string10 string9];
    end
    if i == 3
        string = [string2 string4 string5 string6 string9];
        string_extra = [string1 string4 string5 string6 string10 string9];
    end
    if i == 4
        string = [string2 string4 string5 string6 string7 string9];
        string_extra = [string1 string4 string5 string6 string7 string10 string9];
    end
    if i == 5
        string = [string2 string4 string5 string6 string7 string8 string9];
        string_extra = [string1 string4 string5 string6 string7 string8 string10 string9];
    end
    load(string);
%     if i == 1
%         semilogy(SNR_Values, SNR_ErrorRates','LineWidth',1.5);
%     end
%     if i > 1
%         snrval = SNR_Values;
%         snrerr = SNR_ErrorRates';
%         load(string_extra); 
%         semilogy([snrval SNR_Values], [snrerr SNR_ErrorRates'],'LineWidth',1.5);
%     end
    semilogy(SNR_Values, SNR_ErrorRates','LineWidth',1.5);
end
title('BER for Rayleigh Channel with 1x1 System');
xlabel('SNR Values');
ylabel('Bit Error Rate (BER)');
legend('Ideal Simulation','NonIdeal: AGC','NonIdeal: AGC, BBD','NonIdeal: AGC, BBD, CE','NonIdeal: AGC, BBD, CE, FE','Location','best');
