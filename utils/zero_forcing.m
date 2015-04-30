%% one tap zero-forcing equalization (assuming perfect channel knowledge)
function out = zero_forcing(in,estimateChannel)
    out = zeros(size(in));
    for i = 1:size(in,2)
        out(:,i) = in(:,i)./estimateChannel;
    end
 
