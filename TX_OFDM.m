function  [ofdm, dataMod]  =  TX_OFDM(dataIn, M, N, usedN, CP  )
k  =  log2(M);
unusedN = N-usedN;

dataInMatrix  =  reshape(dataIn, length(dataIn)/k, k);	
dataSymbolsIn  =  bi2de(dataInMatrix);	

dataMod  =  qammod(dataSymbolsIn, M, 'gray');
dataModUsedN  =  reshape(dataMod, usedN, length(dataMod)/usedN);
dataModN  =  zeros(N, size(dataModUsedN, 2));

ofdm  =  zeros(size(dataModN, 2)*(N+CP), 1); 
ofdmSymbol  =  zeros(N, size(dataModN, 2)); 
ofdmSymbolCP  =  zeros(N+CP, size(dataModN, 2));

for  j=1:size(dataModN, 2)
dataModN(:, j) = vertcat(zeros(unusedN/2, 1),dataModUsedN(:, j),zeros(unusedN/2, 1)); 
ofdmSymbol(:, j)  =  ifft(dataModN(:, j), N);
ofdmSymbolCP(:, j)  =  vertcat(ofdmSymbol(N-CP+1:N, j), ofdmSymbol(:, j)); 
ofdm((j-1)*size(ofdmSymbolCP)+1:j*size(ofdmSymbolCP))  =  ofdmSymbolCP(:, j);
end
end
