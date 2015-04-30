function [out,remBits] = bit_power_loading(in, bitAlloc, powerAlloc, modulator)

    nSubcar             = length(bitAlloc);                         % number of Subcarrier
    totalBits           = length(in);                               % total bits need to modulate
    nBitsPerOfdmSymbol  = sum(bitAlloc);                            % number of bits per ofdm symbol
    nOfdmSymbol         = floor(totalBits/nBitsPerOfdmSymbol);      % number of ofdm symbol, floor rounding
    
    
    remBits             = rem(totalBits, nBitsPerOfdmSymbol);       % rem bits not modulated
    
    in = reshape(in(1:nOfdmSymbol*nBitsPerOfdmSymbol),nBitsPerOfdmSymbol,nOfdmSymbol); 
    out= zeros(nSubcar,nOfdmSymbol);

    bitsIndex = 0;
    for i = 1:length(bitAlloc)
        if bitAlloc(i) ~= 0    % if the subcarrier is loading with bits            
            modOrder = 2^bitAlloc(i);
            % for input of QAM modulator, input must be a column vector
            out(i,:) = sqrt(powerAlloc(i))*modulator(reshape(in(bitsIndex+1:bitsIndex+bitAlloc(i),:),bitAlloc(i)*nOfdmSymbol,1), modOrder);
            bitsIndex = bitsIndex + bitAlloc(i);
        end
        
    end
    out = reshape(out,nSubcar*nOfdmSymbol,1);  %% note that the function of reshape is taken columwise, so the symbol order are not corrupted

end