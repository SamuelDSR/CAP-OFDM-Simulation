function constellation = plot_qam_constellation(modOrder)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%% Plot constellation
handle = comm.RectangularQAMModulator(modOrder,'BitInput',true,'SymbolMapping','Gray');
handle.NormalizationMethod = 'Minimum distance between symbols';
bitsperSymbol = log2(modOrder);
int_in = 0:modOrder-1;
bits_in = de2bi(int_in);
bits_in = reshape(bits_in',bitsperSymbol*modOrder,1);
constellation = step(handle,bits_in);
figure
scatterplot(constellation);
if bitsperSymbol <= 5
    text(real(constellation), imag(constellation)+0.1, dec2bin(int_in))
end
title(strcat(num2str(modOrder),'-QAM, Gray Symbol Mapping'))
axis([-modOrder/2 modOrder/2 -modOrder/2 modOrder/2])

end

