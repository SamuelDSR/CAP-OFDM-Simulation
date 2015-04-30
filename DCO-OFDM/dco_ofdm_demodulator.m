function out = dco_ofdm_demodulator(in, demodParams)
nSubcar             = demodParams(1);
cpSize              = demodParams(2);

fftSize     = nSubcar*2+2; % FFT size (this relation is only true in ACO-OFDM)
blkSize     = fftSize+cpSize;  % OFDM symbol length


if mod(length(in),blkSize) == 0   
    nOfdmSymbol = length(in)/blkSize;    
    in          = reshape(in,blkSize,nOfdmSymbol);
    inFFT       = 1/sqrt(fftSize)*fft(in(cpSize+1:end,:));     % cycle prefix remove and conversion to frequency domain
    out         = inFFT(2:nSubcar+1,:);                        % Pilot remove
    out         = reshape(out,nSubcar*nOfdmSymbol,1);
else
    error('DCO-OFDM wrong input size for demodulation');
end


%% Order Version
% totalChannelCoeff   = cell2mat(demodParams(1));
% cpSize              = cell2mat(demodParams(2));
% nSubcar             = cell2mat(demodParams(3));
% 
% fftSize = (nSubcar+1)*2;    % FFT size (this relation is only true in ACO-OFDM)
% blkSize = fftSize+cpSize;   % OFDM symbol length including Cycle Prefix
% 
% if mod(length(in),blkSize) == 0  
%     
%     nOfdmSymbol = length(in)/blkSize;
%     
%     out = zeros(nOfdmSymbol*nSubcar,1);  %pre-allocate memeory
%     
%     for i = 1:nOfdmSymbol
%   
%         inBlk = in(1+(i-1)*blkSize:i*blkSize);
%         
%         % cp removal
%         inCpRemove = inBlk(cpSize+1:end);
%               
%         % normalized fft
%         inFFT = 1/sqrt(fftSize)*fft(inCpRemove);
%         
%         %one tap zero-forcing equalization (assuming perfect channel knowledge)
%         totalChannelFR      = fft(totalChannelCoeff,fftSize);
%         estimateChannelFR   = totalChannelFR(2:nSubcar+1);
%         inFFT(2:nSubcar+1)  = inFFT(2:nSubcar+1)./estimateChannelFR;
%         
%         % pilot removal for DCO-OFDM
%         inPilotRemove = inFFT(2:nSubcar+1);
%         
%         out(1+(i-1)*nSubcar:i*nSubcar) = inPilotRemove;       
%     end
%     
% else
%     error('DCO-OFDM wrong input size for demodulation');
% end