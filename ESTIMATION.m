function [channelCorrection] = ESTIMATION(H, nSymbEst, EbN0_dB, k, N, usedN, CP)
unusedN = N - usedN;
%TX
dataModUsedN = ones(usedN, nSymbEst);
dataModN = zeros(N, nSymbEst);

ofdm = zeros(size(dataModN, 2)*N, 1);
ofdmSymbol = zeros(N, size(dataModN, 2));
ofdmSymbolCP = zeros(N+CP, size(dataModN, 2));

for j=1:size(dataModN, 2)
    dataModN(:, j) = vertcat(zeros(unusedN/2, 1), dataModUsedN(:, j), zeros(unusedN/2, 1));
    ofdmSymbol(:, j) = ifft(dataModN(:, j), N);
    ofdmSymbolCP(:, j) = vertcat(ofdmSymbol(N-CP+1:N, j), ofdmSymbol(:, j));
    ofdm((j-1)*length(ofdmSymbolCP)+1:j*length(ofdmSymbolCP)) = ofdmSymbolCP(:, j);
end
% channel
ofdmChannel = filter(H, 1, ofdm);    
%AWGN
EbN0 = 10^(EbN0_dB/10);
snr = (N/(N+CP))*(usedN/N)*EbN0*k;
snrdB = 10*log10(snr);
ofdmAWGN = awgn(ofdmChannel, snrdB, 'measured');
ofdmSymbolCPRx = reshape(ofdmAWGN, N+CP, length(ofdmAWGN)/(N+CP));
ofdmSymbolRx = ofdmSymbolCPRx(CP+1:end, :);
dataModNRx = zeros(N, size(ofdmSymbolRx, 2));
channelCorrectionMatrix = zeros(N, size(ofdmSymbolRx, 2));
for j=1:size(ofdmSymbolRx, 2)
    dataModNRx(:, j) = fft(ofdmSymbolRx(:, j), N);
    channelCorrectionMatrix(:, j) = dataModN(:, j)./dataModNRx(:, j);
end
channelCorrection = mean(channelCorrectionMatrix, 2);

% channel estimation figure
% figure
% subplot(4, 1, 1)
% stem(real(dataModN(:, 1)))
% title('test data sent')
% subplot(4,1, 2)
% stem(real(dataModNRx(:, 1)))
% title('test data received')
% subplot(4, 1, 3)
% stem(real(channelCorrection))
% title('channel correction')
% subplot(4, 1, 4)
% stem(real(channelCorrection.*dataModNRx(:, 1)))
% title('test data received x channel correction')
end

