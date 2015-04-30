function Pb = theory_ber_rectqam(order,snr)
    if rem(log2(order),2) == 0  %Square QAM
        g = 3/2/(order-1);
    else %Rectangular QAM
        g = 6/(5*order-4);
    end
    
    Pb = 0.2*exp(-g*snr);
 