function  [ofdmChannel]  =  CHANNEL(ofdm, H)
ofdmChannel  =  filter(H, 1, ofdm);	                
end

