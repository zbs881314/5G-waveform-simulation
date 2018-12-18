clear all
close all
% input data
M = 64;         %size of signal constellation
k = log2(M);      % number of bits per symbol
N = 1024;        % number of total carriers
usedN = 600;       %number of data carriers
unusedN = N - usedN;       %number of guard carriers

nSymbOFDM = 1000;    % number of OFDM symbols input
n = usedN*k*nSymbOFDM;      % number of bits

CP = N/8;         
ZT = N/8;        

nSymbEst = 2;      

t = 0:1:120;       
BW = 20e6;   
Ts = (1/BW)*1e9;   


EbN0_dB = 0:1:20;

% Intialization
HdB = -inf.*ones(1, length(t));      
SER = ones(1, length(EbN0_dB)).*inf;             
BER = ones(1, length(EbN0_dB)).*inf;                
BER_ZT = ones(1, length(EbN0_dB)).*inf;         
BER_QAM = ones(1, length(EbN0_dB)).*inf;       

% single-path channel
HdB(1) = 0;

% % [EPA] extended pedestrian A model
% HdB(1) = 0;
% HdB(ceil(51/Ts)) = -1;
% HdB(ceil(71/Ts)) = -2;
% HdB(ceil(91/Ts)) = -3;
% HdB(ceil(111/Ts)) = -8;
% HdB(ceil(191/Ts)) = -17.2;
% HdB(ceil(411/Ts)) = -20.8;

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


% input bits
dataIn = randi([0 1], n, 1);             

% TX
[ofdm, dataMod] = TX_OFDM(dataIn, M, N, usedN, CP);

% NOISE
for z=1:length(EbN0_dB)
    [channelCorrection] = ESTIMATION(H, nSymbEst, EbN0_dB(z), k, N, usedN, CP);
    [ofdmChannel] = CHANNEL(ofdm, H);
    [ofdmAWGN] = AWGN(EbN0_dB(z), ofdmChannel, k, N, usedN, CP);
    [dataInRx, dataModRx] = RX(ofdmAWGN, M, N, usedN, CP, channelCorrection);
    dataSymbolsIn = qamdemod(dataMod, M, 'gray');
    dataSymbolsInRx = qamdemod(dataModRx, M, 'gray');
    [SER(z)] = sum(dataSymbolsIn ~= dataSymbolsInRx)./length(dataSymbolsIn);
    [~, BER(z)] = biterr(dataIn, dataInRx);
end


% TX
[ofdmZT, dataModZT] = TX_ZEROTAIL(dataIn, M, N, usedN, ZT);

% noise
for z=1:length(EbN0_dB)
    [channelCorrectionZT] = ESTIMATION_ZEROTAIL(H, nSymbEst, EbN0_dB(z), k, N, usedN, ZT);
    [ofdmChannelZT] = CHANNEL(ofdmZT, H);
    [ofdmAWGNZT] = AWGN(EbN0_dB(z), ofdmChannelZT, k, N, usedN, ZT);
    [dataInRxZT, dataModRxZT] = RX_ZEROTAIL(ofdmAWGNZT, M, N, usedN, ZT, channelCorrectionZT);
    dataSymbolsZT = qamdemod(dataModZT, M, 'gray');
    dataSymbolsInRxZT = qamdemod(dataModRxZT, M, 'gray');
    [SER_ZT(z)] = sum(dataSymbolsZT ~= dataSymbolsInRxZT)./length(dataSymbolsZT);
    [~, BER_ZT(z)] = biterr(dataIn, dataInRxZT);
end


% QAM
dataInMatrix = reshape(dataIn, length(dataIn)/k, k);
dataSymbolsIn = bi2de(dataInMatrix);
dataMod = qammod(dataSymbolsIn, M, 'gray');

% NOISE
for z=1:length(EbN0_dB)
    snrdB = EbN0_dB(z) + 10*log10(k);
    dataModNoise = awgn(dataMod, snrdB, 'measured');

    dataSymbolsInRx = qamdemod(dataModNoise, M, 'gray');
    dataInMatrixRx = de2bi(dataSymbolsInRx);
    dataInRx = reshape(dataInMatrixRx, size(dataInMatrixRx, 1)*size(dataInMatrixRx, 2), 1);
    [SER(z)] = sum(dataSymbolsIn ~= dataSymbolsInRx)./length(dataSymbolsIn);
    [~, BER_QAM(z)] = biterr(dataIn, dataInRx);
end
% Theoretical BER curve
EbN0 = 10.^(EbN0_dB/10);
SER_MQAM = 2*erfc(sqrt((3*k*EbN0)/(2*(M-1))));
BER_MQAM = SER_MQAM./k;

figure
semilogy(EbN0_dB, [BER; BER_ZT; BER_QAM; BER_MQAM], 'linewidth', 2)
grid on
legend('simulation OFDM', 'simulation OFDM ZT', 'simulation QAM', 'Theoretical');
xlabel('EbN0 (dB)'); ylabel('bit error rate');
title('BER OFDM')