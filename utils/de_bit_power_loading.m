function out = de_bit_power_loading(in, bitAlloc, powerAlloc, demodulator )

nSubcar             = length(bitAlloc);   % number of nSubcar per ofdm symbol
nBitsPerOfdmSymbol  = sum(bitAlloc);      % number of bits per ofdm symbol
nOfdmSymbol         = length(in)/nSubcar; % number of ofdm symbol

in                  = reshape(in,nSubcar,nOfdmSymbol);
out                 = zeros(nBitsPerOfdmSymbol,nOfdmSymbol);

bitsIndex           = 0;

% De bit power loading, different nSubcar has different constellation size and power allocated
for i = 1:nSubcar
    if bitAlloc(i) ~= 0 % if the subcarrier is loading with bits, demodulated, otherwise do nothing
       modOrder = 2^bitAlloc(i); 
       % in of demodulator must be a column vector, (attention!, complex
       % vector a, a' == conj(transpose)), so, reshape instead of ' to
       % transform a row vector to a column vector
       out(bitsIndex+1:bitsIndex+bitAlloc(i),:) = reshape(demodulator((reshape(in(i,:), nOfdmSymbol,1))/sqrt(powerAlloc(i)), modOrder),bitAlloc(i),nOfdmSymbol);
       bitsIndex = bitsIndex + bitAlloc(i);
    end    
end

out = reshape(out,nBitsPerOfdmSymbol*nOfdmSymbol,1);
