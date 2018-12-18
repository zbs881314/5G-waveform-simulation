function  [ofdm, dataMod]  =  TX_ZEROTAIL(dataIn, M, N, usedN, ZT)
k  =  log2(M);
unusedN = N-usedN;

dataInMatrix  =  reshape(dataIn,length(dataIn)/k, k);	
dataSymbolsIn  =  bi2de(dataInMatrix);	

dataMod  =  qammod(dataSymbolsIn, M, 'gray');
dataModUsedN  =  reshape(dataMod, usedN, length(dataMod)/usedN);
dataModN  =  zeros(N, size(dataModUsedN, 2));

ofdm  =  zeros(size(dataModN, 2)*N, 1); 
ofdmSymbol  =  zeros(N, size(dataModN, 2)); 
ofdmSymbolZT  =  zeros(N+ZT,size(dataModN,2 )); 
zerotail  =  zeros(ZT, size(dataModN, 2));

for  j=1:size(dataModN, 2)
dataModN(:, j) = vertcat(zeros(unusedN/2, 1), dataModUsedN(:, j), zeros(unusedN/2, 1)); 
ofdmSymbol(:, j)  =  ifft(dataModN(:, j), N);
ofdmSymbolZT(:, j)  =  vertcat(ofdmSymbol(:, j),zerotail(:, j)); 
ofdm((j-1)*size(ofdmSymbolZT)+1:j*size(ofdmSymbolZT))  =  ofdmSymbolZT(:, j);
end
end
