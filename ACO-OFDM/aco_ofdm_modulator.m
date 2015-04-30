function [out, blkSize] = aco_ofdm_modulator(in, modulator_param)
% get parameters
nSubcar   = modulator_param(1);
cpSize    = modulator_param(2);
fftSize   = nSubcar*4; % it's in ACO-OFDM

if mod(length(in),nSubcar) == 0 
    
    Nsym = length(in)/nSubcar;
    
    blkSize = fftSize+cpSize;
    
    out = zeros(Nsym*blkSize,1);  %pre-allocate memeory
    
    for i = 1:Nsym
  
        in_temps = in(1+(i-1)*nSubcar:i*nSubcar);
        
        % pilot insertion and hermetian sysmetry to ensure a real data after FFT
        pilot_ins_data = zeros(fftSize,1);
        for j = 1:nSubcar
            pilot_ins_data(j*2) = in_temps(j);
            pilot_ins_data(fftSize - j*2+2) = conj(in_temps(j));
        end

        % fourier transform time doamain data
        IFFT_data =sqrt(fftSize)*ifft(pilot_ins_data);
       
        % remove negative part
        IFFT_data(IFFT_data<0) = 0;
        
        % add cycle prefix 
        out((i-1)*blkSize+1:i*blkSize) = [IFFT_data(end-cpSize+1:end); IFFT_data];
        
    end
else
    error('ACO-OFDM wrong input size for modulation');
end
        
        
        