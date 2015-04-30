function out = deblock_transmission(in,blockSize,cpSize)
    nBlock = length(in)/blockSize;
    in = reshape(in, blockSize,nBlock);
    out = in(cpSize+1:end,:);
    out = reshape(out,(blockSize-cpSize)*nBlock,1);
end

