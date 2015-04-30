%% ACO Demodulation 
% in: input time-domain data, 1-D column array
% out: subcarrier data, dimension (nSubcar,nOfdmSymbol)

function out = aco_ofdm_demodulator(in, demodParams)
nSubcar             = demodParams(1);
cpSize              = demodParams(2);

fftSize     = nSubcar*4; % FFT size (this relation is only true in ACO-OFDM)
blkSize     = fftSize+cpSize;  % OFDM symbol length


if mod(length(in),blkSize) == 0
    
    nOfdmSymbol = length(in)/blkSize;    
    in          = reshape(in,blkSize,nOfdmSymbol);
    inFFT       = 1/sqrt(fftSize)*fft(in(cpSize+1:end,:));     % cycle prefix remove and conversion to frequency domain
    out         = 2*inFFT(2:2:2*nSubcar,:);    % Pilot remove and composenate the 1/2 attenuation in ACO-OFDM
    out         = reshape(out,nSubcar*nOfdmSymbol,1);
else
    error('ACO-OFDM wrong input size for demodulation');
end