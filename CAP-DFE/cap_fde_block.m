% forming cap-FDE block
% using method of zero paddings 

function out = cap_fde_block(in, nSymPerBlock, paddingzeros )
   if rem(length(in),nSymPerBlock) == 0
       
       nBlock = length(in)/nSymPerBlock;
       blockSize = nSymPerBlock + paddingzeros;
       
       out = zeros(blockSize, nBlock);       
       in = reshape(in,nSymPerBlock,nBlock);
       
       % cycle prefix
       %out(end-nSymPerBlock+1:end,:) = in;
       %out(1:paddingzeros,:) = 0;
       
       out(1:nSymPerBlock,:) = in;
       out(nSymPerBlock+1:end,:) = 0;
   else
       error('Wrong Input Size!');
   end  
   
   out = reshape(out,blockSize*nBlock,1);
   
   %for cycle prefix, add tail zeros to ensure the composenate the filter delay of the last symbol
   %out(end+1:end+paddingzeros) = 0;
end

