clear all
close all
% input data
M = 64;         %size of signal constellation
k = log2(M);      % number of bits per symbol
N = 1024;        % number of total carriers
usedN = 600;       %number of data carriers
unusedN = N - usedN;       %number of guard carriers

nSymbOFDM = 100;    % number of OFDM symbols input
n = usedN*k*nSymbOFDM;      % number of bits

CP = N/8;         
ZT = N/8;         

nSymbEst = 2;      % number of channel estimation OFDM symbols

EbN0_dB = 30;    
disp('EbN0 dB for all transmission:');
disp(EbN0_dB);

t = 0:1:120;        
BW = 20e6;   
Ts = (1/BW)*1e9;    

% Intialization
HdB = -inf.*ones(1, length(t));      

% single-path channel
HdB(1) = 0;

% % [EPA] extended pedestrian A model
HdB(1) = 0;
HdB(ceil(51/Ts)) = -1;
HdB(ceil(71/Ts)) = -2;
HdB(ceil(91/Ts)) = -3;
HdB(ceil(111/Ts)) = -8;
HdB(ceil(191/Ts)) = -17.2;
HdB(ceil(411/Ts)) = -20.8;

% % [EVA] extended vehicular A model
% HdB(1) = 0;
% HdB(ceil(51/Ts)) = -1.5;
% HdB(ceil(151/Ts)) = -1.4;
% HdB(ceil(311/Ts)) = -3.6;
% HdB(ceil(371/Ts)) = -0.6;
% HdB(ceil(711/Ts)) = -9.1;
% HdB(ceil(1091/Ts)) = -7;
% HdB(ceil(1731/Ts)) =-12;
% HdB(ceil(2511/Ts)) = -16.9;

% % [EVA] extended typical Urban model
% HdB(1) = -1;
% HdB(ceil(51/Ts)) = -1;
% HdB(ceil(121/Ts)) = -1;
% HdB(ceil(201/Ts)) = 0;
% HdB(ceil(231/Ts)) = 0;
% HdB(ceil(501/Ts)) = 0;
% HdB(ceil(1601/Ts)) = -3;
% HdB(ceil(2301/Ts)) =-5;
% HdB(ceil(5001/Ts)) = -7;


H = 10.^(HdB/10);              

figure
stem(t, HdB)
title('Channel')
xlabel('Time Sample')
ylabel('Signal realtion (dB)')

% Theoretical bit error rate
EbN0 = 10.^(EbN0_dB/10);
SER_MQAM = 2*erfc(sqrt((3*k*EbN0)/(2*(M-1))));
BER_MQAM = SER_MQAM./k;

disp('Theoretical BER')
disp(BER_MQAM);

% input bits
dataIn = randi([0 1], n, 1);    


% channel estimation
[channelCorrection] = ESTIMATION(H, nSymbEst, EbN0_dB, k, N, usedN, CP);

% TX
[ofdm, dataMod] = TX_OFDM(dataIn, M, N, usedN, CP);

% CHANNEL
[ofdmChannel] = CHANNEL(ofdm, H);

% NOISE
[ofdmAWGN] = AWGN(EbN0_dB, ofdmChannel, k, N, usedN, CP);

% RX
[dataInRx, dataModRxFixed] = RX(ofdmAWGN, M, N, usedN, CP, channelCorrection);

% BER
[~, BER] = biterr(dataIn, dataInRx);          
disp('BER OFDM CHANNEL FIXED');
disp(BER);

% Constellation
sPlotFig = scatterplot(dataModRxFixed, 1, 0, 'g.');
grid on
hold on
scatterplot(dataMod, 1, 0, 'k*', sPlotFig)
title('OFDM AWGN CHANNEL CONSTELLATION')

% channel estimation
[channelCorrection] = ESTIMATION_ZEROTAIL(H, nSymbEst, EbN0_dB, k, N, usedN, ZT);

% TX
[ofdmZT, dataMod] = TX_ZEROTAIL(dataIn, M, N, usedN, ZT);

% channel
[ofdmChannel] = CHANNEL(ofdmZT, H);

% noise
[ofdmAWGN] = AWGN(EbN0_dB, ofdmChannel, k, N, usedN, ZT);

% RX
[dataInRx, dataModRxFixed] = RX_ZEROTAIL(ofdmAWGN, M, N, usedN, ZT, channelCorrection);

% BER
[~, BER] = biterr(dataIn, dataInRx);
disp('BER OFDM zerotail channel fixed');
disp(BER);

% constellation
sPlotFig = scatterplot(dataModRxFixed, 1, 0, 'g.');
hold on
grid on
scatterplot(dataMod, 1, 0, 'k*', sPlotFig)
title('OFDM zerotail channel fixed constellation')

% QAM
dataInMatrix = reshape(dataIn, length(dataIn)/k, k);
dataSymbolsIn = bi2de(dataInMatrix);
dataMod = qammod(dataSymbolsIn, M, 'gray');

% NOISE
snrdB = EbN0_dB + 10*log10(k);
%AWGN channel
dataModNoise = awgn(dataMod, snrdB, 'measured');

% demodulate with Gray code
dataSymbolsInRx = qamdemod(dataModNoise, M, 'gray');
% convert to binary data
dataInMatrixRx = de2bi(dataSymbolsInRx);
%reshape into a binary vector
dataInRx = reshape(dataInMatrixRx, size(dataInMatrixRx, 1)*size(dataInMatrixRx, 2), 1);
% BER
[~, BER] = biterr(dataIn, dataInRx);
disp('BER QAM');
disp(BER);

% constellation
sPlotFig = scatterplot(dataModNoise, 1, 0, 'g.');
hold on
grid on
scatterplot(dataMod, 1, 0, 'k*', sPlotFig)
title('QAM constellation')