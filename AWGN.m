function [ofdmAWGN] = AWGN(EbN0_dB, ofdm, k, N, usedN, CP)
EbN0 = 10^(EbN0_dB/10);
snr = (N/(N+CP))*(usedN/N)*EbN0*k;
snrdB = 10*log10(snr);
ofdmAWGN = awgn(ofdm, snrdB, 'measured');

% AWGN addition figure
% figure
% plot(real(ofdmAWGN), 'r')
% hold on
% plot(real(ofdm))
% plot(real(ofdmAWGN)-real(ofdm), 'g')
% title('OFDM NOISE')
% legend('OFDM AWGN', 'OFDM', 'NOISE')

end
