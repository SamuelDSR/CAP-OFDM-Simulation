function out = led_shaping(in, led)
%get dc bias
out = in + led.dcBias;

%polynomial distortion
out = led_poly_non_linear_shaping(out,led.polynomial);
 
%clipping 
out(out<led.minCurrent) = led.minCurrent;
out(out>led.maxCurrent) = led.maxCurrent;
