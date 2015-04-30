function out = led_poly_non_linear_shaping(in, polynomial)
    poly_length = length(polynomial);
    out = zeros(size(in));
    for i = 1:poly_length
        out = out+polynomial(i)*in.^(i-1);
    end