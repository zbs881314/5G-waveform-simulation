function  [ofdm, dataMod]  =  TX_DFT(dataIn, M, N, usedN, ZT)
k  =  log2(M);
unusedN = N-usedN;
ZTh = 10;	       
ZTt = ZT-ZTh;	 

dataInMatrix  =  reshape(dataIn, length(dataIn)/k, k);	
dataSymbolsIn  =  bi2de(dataInMatrix);	

dataMod  =  qammod(dataSymbolsIn, M, 'gray');
reshapeUsedN  =  floor(length(dataMod)/(usedN-ZT));	
disp('Discarded  last  QAM  symbols  in  ZT©\DFTs');
disp(length(dataMod)-reshapeUsedN*(usedN-ZT)); 
dataMod = dataMod(1:reshapeUsedN*(usedN-ZT));

dataModUsedN  =  reshape(dataMod, usedN-ZT, reshapeUsedN);	
dataModN  =  zeros  (N, size(dataModUsedN, 2));	

ofdm  =  zeros(size(dataModN, 2)*(N+ZT), 1);
ofdmSymbol  =  zeros(usedN-ZTh-ZTt, size(dataModN, 2)); 
ofdmSymbolZT  =  zeros(usedN, size(dataModN, 2)); 
ofdmSymbolDFT  =  zeros(N, size(dataModN, 2));

for  j=1:size(dataModN, 2)
ofdmSymbol(:, j)  =  ifft(dataModUsedN(:, j), usedN-ZT);
ofdmSymbolZT(:, j)  =  vertcat(zeros(ZTh,1), ofdmSymbol(:, j), zeros(ZTt, 1));
dataModUsedNZT  =  fft(ofdmSymbolZT, usedN);
dataModN(:,j) = vertcat(zeros(unusedN/2, 1), dataModUsedNZT(:, j), zeros(unusedN/2, 1));
ofdmSymbolDFT(:, j)  =  ifft(dataModN(:, j), N); 
ofdm((j-1)*size(ofdmSymbolDFT)+1:j*size(ofdmSymbolDFT))  =  ofdmSymbolDFT(:, j);
end
end

