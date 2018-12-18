function  [channelCorrection]  =  ESTIMATION_ZEROTAIL(H,  nSymbEst,  EbN0_dB,  k,  N, usedN,  ZT)

unusedN = N - usedN;

%TX
dataModUsedN  =  ones(usedN, nSymbEst);	    
dataModN  =  zeros(N, nSymbEst);	              

ofdm  =  zeros(size(dataModN, 2)*N, 1); 
ofdmSymbol  =  zeros(N, size(dataModN, 2)); 
ofdmSymbolZT  =  zeros(N+ZT, size(dataModN, 2)); 
zerotail  =  zeros(ZT, size(dataModN, 2));

for  j=1:size(dataModN, 2)
dataModN(:, j) = vertcat(zeros(unusedN/2, 1), dataModUsedN(:, j), zeros(unusedN/2, 1));
ofdmSymbol(:, j)  =  ifft(dataModN(:, j), N);
ofdmSymbolZT(:, j)  =  vertcat(ofdmSymbol(:, j), zerotail(:, j));
ofdm((j-1)*length(ofdmSymbolZT)+1:j*length(ofdmSymbolZT))  =  ofdmSymbolZT(:, j);
end

%CHANNEL
ofdmChannel  =  filter(H, 1, ofdm);	  

%AWGN
EbN0 = 10^(EbN0_dB/10);
snr = (N/(N+ZT))*(usedN/N)*EbN0*k;
snrdB  =  10*log10(snr);	
ofdmAWGN  =  awgn(ofdmChannel, snrdB, 'measured');	
%RX
ofdmSymbolZTRx  =  reshape(ofdmAWGN, N+ZT, length(ofdmAWGN)/(N+ZT)); 
dataModNRx  =  zeros(N, size(ofdmSymbolZTRx, 2)); 
channelCorrectionMatrix  =  zeros(N, size(ofdmSymbolZTRx, 2));

for  j=1:size(ofdmSymbolZTRx, 2)
ofdmSymbolZTRx(1:ZT, j)  =  ofdmSymbolZTRx(1:ZT, j)  +  ofdmSymbolZTRx(N+1:end, j); 
ofdmSymbolRx  =  ofdmSymbolZTRx(1:N, :);
dataModNRx(:, j)  =  fft(ofdmSymbolRx(:, j), N);
channelCorrectionMatrix(:, j)=dataModN(:, j)./dataModNRx(:, j);
end

channelCorrection  =  mean(channelCorrectionMatrix, 2);
%  Channel  estimation  figure
% figure
% subplot(4,1,1)
%  stem(real(dataModN(:,1)))
%  title('Test  Data  Sent')
% subplot(4,1,2)
%  stem(real(dataModNRx(:,1)))
%  title('Test  Data  Received')
% subplot(4,1,3)
%  stem(real(channelCorrection))
%  title('Channel  Correction')
% subplot(4,1,4)
%  stem(real(channelCorrection.*dataModNRx(:,1)))
%  title('Test  Data  Received  X  Channel  Correction')

end


