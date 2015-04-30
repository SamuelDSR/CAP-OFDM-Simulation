function out = block_transmission(in, nSymPerBlock, cpSize )
   if rem(length(in),nSymPerBlock) == 0
       nBlock = length(in)/nSymPerBlock;
       blockSize = nSymPerBlock + cpSize;
       
       out = zeros(blockSize, nBlock);       
       in = reshape(in,nSymPerBlock,nBlock);
       
       out(end-nSymPerBlock+1:end,:) = in;
       out(1:cpSize,:) = out(end-cpSize+1:end,:);
   else
       error('Wrong Input Size!');
   end  
   
   out = reshape(out,blockSize*nBlock,1);
end

