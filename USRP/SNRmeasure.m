N = 10^(-7.86);
SNdb = [19 8 20 10 30 18 18 25 25 20 32 29 28 28 25 26 26 20 26 31]-78.6;
SN = 10.^(SNdb/10);
S = SN-N;
SNR = 10*log10(S/N)