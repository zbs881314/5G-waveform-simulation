function  [dataInRx, dataModRxFixed]  =  RX_ZEROTAIL(ofdm, M, N, usedN, ZT, channelCorrection)

unusedN = N-usedN;

ofdmSymbolZTRx  =  reshape(ofdm, N+ZT, length(ofdm)/(N+ZT));

dataModNRx  =  zeros(N, size(ofdmSymbolZTRx, 2)); 
dataModUsedNRx  =  zeros(usedN, size(ofdmSymbolZTRx, 2)); 
dataModNRxFixed  =  zeros(N, size(ofdmSymbolZTRx, 2));

for  j=1:size(ofdmSymbolZTRx, 2) 
ofdmSymbolZTRx(1:ZT, j)  =  ofdmSymbolZTRx(1:ZT, j)  +  ofdmSymbolZTRx(N+1:end, j); 
ofdmSymbolRx  =  ofdmSymbolZTRx(1:N, j);
dataModNRx(:, j)  =  fft(ofdmSymbolRx, N);
dataModNRxFixed(:, j)  =  dataModNRx(:, j).*channelCorrection;
dataModUsedNRx(:, j)  =  dataModNRxFixed(unusedN/2+1:N-unusedN/2, j);
end

%  Info  of  the  OFDM  carriers  to  a  chain  of  QAM  simbols  to  get  data  vector 
dataModRxFixed  =  reshape(dataModUsedNRx, size(dataModUsedNRx, 1)*size(dataModUsedNRx, 2), 1); 
dataSymbolsInRx  =  qamdemod(dataModRxFixed, M, 'gray');
dataInMatrixRx  =  de2bi(dataSymbolsInRx); dataInRx  =  reshape(dataInMatrixRx, size(dataInMatrixRx, 1)*size(dataInMatrixRx, 2), 1);

end

