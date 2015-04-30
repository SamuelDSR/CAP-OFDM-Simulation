function out = dotprod(in1, in2)
    [m1,n1] = size(in1);
    [m2,n2] = size(in2);
    if (m1==m2 && n1==n2)
        out = in1.*in2;
    elseif (m1 == n2 && m2 == n1)
        out = in1.*(in2)';
    else
        error('dotprod, array size not match');
    end
    
    if m1<n1
        out = reshape(out,n1,m1);
    end
     
           