% in: 1-D array or 2-D array
% blocksize:  number of symbols per block
% cpsize:     number of cycle prefix
% out: 1-D array

function out = cap_fde_deblock(in,blockSize,cpSize)

    nBlock = length(in)/blockSize;
    in = reshape(in, blockSize,nBlock);
    
    out = in(cpSize+1:end,:);
    out = reshape(out,(blockSize-cpSize)*nBlock,1);
end

