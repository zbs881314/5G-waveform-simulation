clear all
close all
% input data
M = 64;         %size of signal constellation
k = log2(M);      % number of bits per symbol
N = 2048;        % number of total carriers
usedN = 1200;       %number of data carriers
unusedN = N - usedN;       %number of guard carriers

nSymbOFDM = 100;    % number of OFDM symbols input
n = usedN*k*nSymbOFDM;      % number of bits

CP = N/8;         % cyclic prefix length samples
ZT = N/8;         % zero tail length samples

dataIn = randi([0 1], n, 1);

% TX
[ofdmZTDFT, ~] = TX_OFDM_ZEROTAIL_DFT(dataIn, M, N, usedN, ZT);

% TX
[ofdmZT, ~] = TX_OFDM_ZEROTAIL(dataIn, M, N, usedN, ZT);

% TX
[ofdm, ~] = TX_OFDM(dataIn, M, N, usedN, CP);

figure
subplot(3, 1, 1)
plot(real(ofdm(1:N+300)))
title('OFDM')
subplot(3, 1, 2)
plot(real(ofdmZT(1:N+300)), 'r')
title('OFDM ZT')
subplot(3, 1, 3)
plot(real(ofdmZTDFT(1:N+300)), 'g')
title('OFDM ZT DFTs')

figure
[pxx, f] = periodogram(ofdmZT);
plot(f/pi, 10*log10(pxx), 'r')
hold on
[pxx, f] = periodogram(ofdm);
plot(f/pi, 10*log10(pxx))
hold on
[pxx, f] = periodogram(ofdmZTDFT);
plot(f/pi, 10*log10(pxx), 'g')
ylabel('power/frequency (dB/rad/sample)')
xlabel('normalized frequency (xpi rad/sample)')
title('OFDM ZEROTAIL vs DFTs')
legend('OFDM ZT', 'OFDM', 'OFDM ZT DFTs')




