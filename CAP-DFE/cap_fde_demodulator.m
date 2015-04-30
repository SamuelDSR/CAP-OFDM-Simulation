  function out = cap_fde_demodulator(in, nSymbolPerBlock, FilterSpan, UpSamplingFactor, complexMatchFilter, channelStateInformation);
      
      fftSize = length(channelStateInformation);
      blkSize = fftSize;
      
      totFilter  = dotprod(fft(complexMatchFilter', fftSize), 1./channelStateInformation);
      
      in = reshape(in,blkSize,length(in)/blkSize);
      [nSamples, nBlock] = size(in);
      filterDelay = UpSamplingFactor*FilterSpan;
      
      out = zeros(nSymbolPerBlock,nBlock);
                
      for i = 1:nBlock
          % match filtering and equalization in frequency domain
          tmpf  = dotprod(fft(in(:,i),fftSize),totFilter);
          
          % downsampling in frequency domain
%           tmpf  = reshape(tmpf,length(tmpf)/UpSamplingFactor,UpSamplingFactor);        
%           tmpdf = zeros(length(tmpf)/UpSamplingFactor,1);
%           for j = 1:UpSamplingFactor
%               tmpdf = tmpdf+ tmpf(:,j);
%           end

          % downsampling in time domain
          tmpt = ifft(tmpf,fftSize);
          for j = 1:nSymbolPerBlock             
            out(j,i) = tmpt(1+filterDelay+(j-1)*UpSamplingFactor);
          end
      end
      
      out = reshape(out,nSymbolPerBlock*nBlock,1);
end
          
          
      