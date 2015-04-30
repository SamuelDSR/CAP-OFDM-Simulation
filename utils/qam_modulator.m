function out = qam_modulator(in, modOrder)
%% generator QAM modulator
% handle = modem.qammod(order);
% handle.InputType = 'Bit';
% handle.SymbolOrder = 'Gray';
% handle.NormalizationMethod = 'Average power';

%% rectangular QAM modulator
handle = comm.RectangularQAMModulator(modOrder,'BitInput',true,'SymbolMapping','Gray');
handle.NormalizationMethod = 'Average power';

out = step(handle,in);
