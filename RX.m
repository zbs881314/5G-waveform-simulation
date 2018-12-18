function  [dataInRx, dataModRxFixed]  =  RX(ofdm, M, N, usedN, CP, channelCorrection) 
unusedN = N-usedN;
ofdmSymbolCPRx  =  reshape(ofdm, N+CP, length(ofdm)/(N+CP)); 
ofdmSymbolRx  =  ofdmSymbolCPRx(CP+1:end, :);

dataModNRx  =  zeros(N, size(ofdmSymbolRx, 2)); 
dataModUsedNRx  =  zeros(usedN, size(ofdmSymbolRx, 2)); 
dataModNRxFixed  =  zeros(N, size(ofdmSymbolRx, 2));

for  j=1:size(ofdmSymbolRx, 2)
dataModNRx(:, j)  =  fft(ofdmSymbolRx(:, j), N);
dataModNRxFixed(:, j)  =  dataModNRx(:,j ).*channelCorrection;
dataModUsedNRx(:, j)  =  dataModNRxFixed(unusedN/2+1:N-unusedN/2, j);
end

%  Info  of  the  OFDM  carriers  to  a  chain  of  QAM  simbols  to  get  data  vector 
dataModRxFixed  =  reshape(dataModUsedNRx, size(dataModUsedNRx, 1)*size(dataModUsedNRx, 2), 1); 
dataSymbolsInRx  =  qamdemod(dataModRxFixed, M,'gray');
dataInMatrixRx  =  de2bi(dataSymbolsInRx); 
dataInRx  =  reshape(dataInMatrixRx, size(dataInMatrixRx, 1)*size(dataInMatrixRx, 2), 1);
end
