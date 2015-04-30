function out = qam_demodulator(in, modOrder)
%% generator QAM Demodulator
% handle = modem.qamdemod(order);
% handle.OutputType = 'Bit';
% handle.SymbolOrder = 'Gray';
% handle.NormalizationMethod = 'Average power';

%% rectangular QAM Demodulator
handle = comm.RectangularQAMDemodulator(modOrder,'BitOutput',true,'SymbolMapping','Gray');
handle.NormalizationMethod = 'Average power';

out = step(handle,in);
